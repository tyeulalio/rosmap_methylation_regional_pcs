library(RNOmni)
library(DMRcate)
library(minfi)
library(GenomicRanges)
library(limma)
library(ggplot2)
library(tidyverse)

# perform DM and DMR analysis together
# create mixed data frame of true positive and negative genes

savedir <- paste0("../../output/05_combined_dmrs_dmp/")
dir.create(savedir)

# use the colors defined in the color panels R script
source("../color_panels.R")


get_parameters <- function(parameter_choices, num_regions, repeat_genes=FALSE){
    # get a set of parameters
    parameter_rows <- sample(1:nrow(parameter_choices),
                                size=num_regions,
                                replace=TRUE
    ) 
    parameters <- do.call(cbind, 
                             parameter_choices[parameter_rows,,drop=F]) %>%
        as.data.frame()
    head(parameters)
    dim(parameters)

    # choose the gene that we will use from each 
    # ensure that we never choose the same gene multiple times
    # for one set of parameters
    params <- parameters %>%
        group_by(percent_meth_difference,
                 percent_sites_dm,
                 num_sites) %>%
        mutate(gene_num = sample(1:1000, n(), replace=repeat_genes)) %>%
        arrange(percent_meth_difference,
                 percent_sites_dm,
                 num_sites) 
    head(params)
    dim(params)
    params
}


get_runname <- function(params){
    percent_meth_difference <- params[['percent_meth_difference']] %>% as.numeric()
    percent_sites_dm <- params[['percent_sites_dm']] %>% as.numeric()
    num_sites <- params[['num_sites']] %>% as.numeric()
    num_samples <- params[['num_samples']] %>% as.numeric()


    dmr_length <- num_sites / 2
    N <- 1000

    # name of output files are formatted like this
    (runname <- paste0(
                        "_numSites", num_sites,
                        "_numSamples", num_samples,
                        "_pct_meth_diff", percent_meth_difference,
                        "_dmr_length", dmr_length,
                        "_pct_sites_dm", percent_sites_dm,
                        "_N", N))
    runname
}


load_methylation <- function(ordered_parameters, num_samples){
    # load the methylation data that matches the parameters
    head(ordered_parameters)

    # need to load avgs, rpcs, cpgs for each set of parameters
    # then select the genes from that set 
    unique_params <- ordered_parameters %>%
        group_by(percent_meth_difference, percent_sites_dm, num_sites) %>%
        summarize(gene_num=list(gene_num)) %>%
        unique() %>%
        mutate(num_samples = num_samples)
    head(unique_params)
    dim(unique_params)

    params <- unique_params[1,]
    i <- 1 
    load_parameters <- function(i){
        params <- unique_params[i,]
        (runname <- get_runname(params))
        #print(paste(runname, "run", i, "out of", nrow(unique_params)))

        datadir <- "../../output/00_simulated_data/"
        (datafile <- paste0(datadir,
                           "simulated_data",
                           runname, 
                           ".rds"
        ))
        simdata <- readRDS(datafile)
        names(simdata)

        percent_meth_difference <- params[['percent_meth_difference']] %>% as.numeric()
        percent_sites_dm <- params[['percent_sites_dm']] %>% as.numeric()
        num_sites <- params[['num_sites']] %>% as.numeric()
        num_samples <- params[['num_samples']] %>% as.numeric()
        
        # get the genes for this parameter
        gene_nums <- params[['gene_num']][[1]] %>%
            sort() %>%
            unique()
        head(gene_nums)
        #print(gene_nums)
        params
        avgs <- simdata[['avgs']][gene_nums,,drop=FALSE] %>%
            as.data.frame() %>%
            mutate(percent_meth_difference=percent_meth_difference,
                   percent_sites_dm=percent_sites_dm,
                   num_sites=num_sites,
                   num_samples=num_samples
            )
        head(avgs)

        # get rpcs
        rpcs <- simdata[['rpcs']] %>%
            na.omit() %>%
            rownames_to_column('region_pc') %>%
            separate(region_pc, c('region', 'pc'), sep='-', remove=FALSE) %>%
            mutate(region = str_remove(region, 'region'))
        rpcs <- rpcs[match(gene_nums, rpcs$region),] %>%
            select(-region, -pc) %>%
            as.data.frame() %>%
            mutate(percent_meth_difference=percent_meth_difference,
                   percent_sites_dm=percent_sites_dm,
                   num_sites=num_sites,
                   num_samples=num_samples
            )
        head(rpcs)[,1:10]
        colnames(rpcs)

        # get cpgs
        names(simdata)
        raw_data <- simdata[['raw_data']]
        head(raw_data[[1]])
        i <- 10
        get_cpgs <- function(i){
            if (all(is.na(raw_data[[i]]))) return(NA)

            cpgs <- raw_data[[i]][[1]] %>%
                as.data.frame() %>%
                mutate(region = i) %>%
                select(region, cpg_pos=chrpos, starts_with('PM'))
            head(cpgs)[,1:10]
            ncol(cpgs)

            # fix column names 
            samples <- data.frame(samples=colnames(cpgs)[3:ncol(cpgs)]) %>%
                separate(samples, c('pm', 'sampnum'), sep='\\.') %>%
                mutate(sampnum = as.numeric(sampnum)) %>%
                mutate(pheno = ifelse(sampnum <= (num_samples/2), 'case', 'control')) %>%
                mutate(sample_name = paste0(pheno, "_", sampnum))
            samples

            colnames(cpgs)[3:ncol(cpgs)] <- samples$sample_name

            cpgs <- cpgs 
            head(cpgs)[,1:10]
            cpgs
        }
        cpgs <- lapply(gene_nums, get_cpgs) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            mutate(percent_meth_difference=percent_meth_difference,
                   percent_sites_dm=percent_sites_dm,
                   num_sites=num_sites,
                   num_samples=num_samples
            )
        head(cpgs)
        tail(cpgs)
        dim(cpgs)
        length(unique(cpgs$region))
        
        list(avgs=avgs, rpcs=rpcs, cpgs=cpgs, params=params)
    }
    meth <- lapply(1:nrow(unique_params), load_parameters)

    head(ordered_parameters)
    head(unique_params)

    # combine avgs together, remove NAs
    unique_avgs <- lapply(meth, function(x) x[['avgs']]) %>%
        do.call(rbind, .) %>%
        rownames_to_column('region') %>%
        na.omit() %>%
        mutate(gene_num = as.numeric(str_remove(region, 'region')))
    dim(unique_avgs)
    head(unique_avgs)
    head(ordered_parameters)
    avgs <- ordered_parameters %>%
        inner_join(unique_avgs)
    head(avgs)[,1:10]
    head(avgs)
    dim(avgs)

    head(ordered_parameters)
    unique_rpcs <- lapply(meth, function(x) x[['rpcs']]) %>%
        do.call(rbind, .) %>%
        separate(region_pc, c('region', 'pc'), sep='-') %>%
        mutate(gene_num = as.numeric(str_remove(region, 'region'))) %>%
        na.omit()
    rpcs <- ordered_parameters %>%
        inner_join(unique_rpcs)
    head(rpcs)[,1:10]
    head(rpcs)
    dim(rpcs)

    unique_cpgs <- lapply(meth, function(x) x[['cpgs']]) %>%
        do.call(rbind, .) %>%
        rename(gene_num=region) %>%
        na.omit()
    head(unique_cpgs)[,1:10]
    head(unique_cpgs)
    head(ordered_parameters)
    unique_cpgs[14521,]
    cpgs <- ordered_parameters %>%
        left_join(unique_cpgs, relationship='many-to-many')
    head(cpgs)[,1:10]
    head(cpgs)
    dim(cpgs)

    meth <- list(avgs=avgs, rpcs=rpcs, cpgs=cpgs)
}

normalize_meth <- function(meth){
    # normalize the methylation data
    names(meth)

    # normalize avgs
    head(meth$avgs)
    colnames(meth$avgs)

    df <- meth$cpgs
    get_normalized_meth <- function(df){
        # get the methylation data only
        full_df <- na.omit(df)

        df_meth <- full_df[,grepl('case|control', colnames(full_df))] 
        df_params <- full_df[,!grepl('case|control', colnames(full_df))] 
        colnames(df_meth)
        colnames(df_params)

        # remove zero variance cpgs
        # this is a standard practice
        keep_rows <- apply(df_meth, 1, var, na.rm=TRUE) != 0
        var_df <- df_meth[keep_rows,]
        normalized_df <- apply(var_df, 1, RankNorm) %>%
            t() %>%
            as.data.frame()
        head(normalized_df)

        # connect the parameters back
        keep_params <- df_params[keep_rows,]
        normalized_df <- cbind(normalized_df, keep_params)
        normalized_df
    }
    normalized_avgs <- get_normalized_meth(meth$avgs)
    normalized_rpcs <- get_normalized_meth(meth$rpcs)
    normalized_cpgs <- get_normalized_meth(meth$cpgs)

    normalized_meth <- list(avgs=normalized_avgs,
                            rpcs=normalized_rpcs,
                            cpgs=normalized_cpgs
    )
    normalized_meth
}

create_methylation_data <- function(total_genes, true_dm_proportion, num_samples){
    # create data frame of methylation data
    # using mixed parameters

    (true_dm_count <- total_genes * true_dm_proportion)
    (no_dm_count <- total_genes - true_dm_count)

    # options for parameters
    percent_meth_difference_range <- c(seq(0.00001, 0.1, 0.01))
    percent_sites_dm_range <- seq(0.25,0.75,0.25)
    num_sites_range=c(20, 50)
    #num_samples_range=c(50,500)
    #num_samples_range=c(50,500,5000)

    # need to fix the number of samples
    # to do DM on the matrix
    num_samples_ranges <- c(num_samples)

    # randomly select the parameters for the true 
    dm_parameter_choices <- expand.grid(percent_meth_difference=percent_meth_difference_range,
                                     percent_sites_dm=percent_sites_dm_range,
                                     num_sites=num_sites_range
    )
    head(dm_parameter_choices)
    dim(dm_parameter_choices)
    dm_parameters <- get_parameters(dm_parameter_choices, true_dm_count)

    # select the parameter choices for the false
    nodm_parameter_choices <- expand.grid(percent_meth_difference=c(0),
                                     percent_sites_dm=percent_sites_dm_range,
                                     num_sites=num_sites_range
    )
    head(nodm_parameter_choices)
    dim(nodm_parameter_choices)
    nodm_parameters <- get_parameters(nodm_parameter_choices, no_dm_count, repeat_genes=TRUE)
    head(nodm_parameters)
    dim(nodm_parameters)

    # combine the parameters and randomly order them
    # this will be the order of the genes for analysis
    combined_parameters <- rbind(dm_parameters, nodm_parameters)
    head(combined_parameters)
    dim(combined_parameters)
    ordered_parameters <- combined_parameters[sample(1:nrow(combined_parameters), nrow(combined_parameters), replace=FALSE),] 
    head(ordered_parameters)    
    dim(ordered_parameters)    

    # add a chromosome number and number of gene order on chromosome
    ordered_parameters$chr <- sample(1:22, nrow(ordered_parameters), replace=TRUE)
    ordered_parameters <- ordered_parameters %>%
        arrange(chr) %>%
        group_by(chr) %>%
        mutate(gene_order = row_number()) %>%
        ungroup()
    head(ordered_parameters)
    table(ordered_parameters$percent_meth_difference)

    # load the methylation data
    meth <- load_methylation(ordered_parameters, num_samples)
    meth
}

format_meth <- function(normalized_meth){
    # format the meth data for dm analysis
    names(normalized_meth)

    format_gene_level <- function(df){
        head(df)

        num_samples <- df[1,]$num_samples

        # get the samples ordered
        ordered_samples <- data.frame(samples=colnames(df)[grepl('case|control', colnames(df))]) %>%
            separate(samples, c('class', 'sampnum'), sep='_', remove=FALSE) %>%
            mutate(sampnum = as.numeric(sampnum)) %>%
            arrange(sampnum) %>%
            mutate(class = ifelse(sampnum <= (num_samples/2), 'case', 'control')) %>%
            mutate(new_name = paste0(class, "_", sampnum))
        head(ordered_samples)

        df_meth <- df[ordered_samples$samples] 
        colnames(df_meth)
        colnames(df_meth) <- ordered_samples$new_name
        colnames(df_meth)
        head(df_meth)

        df_params <- df[!colnames(df) %in% ordered_samples$samples]
        head(df_params)

        list(meth=df_meth, params=df_params)
    }

    # format gene-level data
    # avgs
    df <- normalized_meth[['avgs']] %>%
        arrange(chr, gene_order) %>%
        mutate(gene = paste0("chr", chr, "_", gene_order)) %>%
        column_to_rownames('gene')
    head(df)
    rownames(df)
    avgs <- format_gene_level(df)

    # rpcs
    df <- normalized_meth[['rpcs']] %>%
        arrange(chr, gene_order, pc) %>%
        mutate(gene = paste0("chr", chr, "_", gene_order, "_", pc)) %>%
        column_to_rownames('gene')
    head(df)
    rownames(df)
    rpcs <- format_gene_level(df)
    head(rpcs)

    df <- normalized_meth[['cpgs']] %>%
        arrange(chr, gene_order, cpg_pos) %>%
        mutate(gene = paste0("chr", chr, "_", gene_order, "_", cpg_pos)) %>%
        column_to_rownames('gene')
    head(df)
    rownames(df)
    cpgs <- format_gene_level(df)
    head(cpgs)
    names(cpgs)
    rownames(cpgs$meth)

    formatted_meth <- list(avgs=avgs, rpcs=rpcs, cpgs=cpgs)
}

#meth <- formatted_meth[['avgs']]
run_dm_analysis <- function(meth){
    # run DM analysis
    head(meth)

    # make sure meth and pheno match
    pheno <- data.frame(samples=colnames(meth)) %>%
        separate(samples, c('class', 'sampnum'), sep='_', remove=FALSE) %>%
        mutate(status = ifelse(class == 'control', 0, 1))
    head(pheno)    

    # run limma dm analysis

    # set up the model + design matrix
    model_content <- paste0("~ status")
    design <- model.matrix(as.formula(model_content), data=pheno)
    
    # run linear model
    lmfit <- lmFit(meth, design)
    
    # ebayes shrinkage
    lm <- eBayes(lmfit)
    
    # grab results
    head(lm)
    results <- topTable(lm, 
                        coef='status',
                        number = Inf,
                        adjust.method='none')
    head(results)
    
    # add p-value adjustment here
    results <- results %>%
        mutate(bh_pval = p.adjust(P.Value, method='BH'),
               bf_pval = p.adjust(P.Value, method='bonferroni'))
    head(results)

    # how many are significant?
    num_sig <- nrow(results[results$bh_pval < 0.05,])
    print(paste("Number of BH sig results:", num_sig))
    
    num_sig <- nrow(results[results$bf_pval < 0.05,])
    print(paste("Number of BF sig results:", num_sig))

    results
}

#dmp_res <- limma_res[['cpgs']]
run_dmr_analysis <- function(dmp_res, fdr_thresh){
    # run dmr analysis with dmrcate
    head(dmp_res)

    meth_annotated <- dmp_res %>%
        mutate(gene_pos = rownames(dmp_res)) %>%
        separate(gene_pos, c('chr', 'region', 'pos'), sep='_') %>%
        select(stat=t, chr, region, pos, diff=logFC, fdr=bh_pval) %>%
        mutate(pos = as.numeric(pos),
               region = as.numeric((region))
        ) %>%
        arrange(chr, region, pos)
    head(meth_annotated)

    # find first region of each chromosom
    head(meth_annotated)
    first_regions <- meth_annotated %>%
        select(chr, region, pos) %>%
        unique() %>%
        group_by(chr, region) %>%
        # find max position for this chromosome and region
        summarize(max_pos = max(pos)) %>%
        mutate(prev_max_pos = lag(max_pos)) %>%
        replace_na(list(prev_max_pos=0)) %>%
        mutate(cumsum_pos = cumsum(prev_max_pos)) %>%
        group_by(chr) %>%
        mutate(first_region = min(region)) %>%
        ungroup()
    first_regions

    head(meth_annotated)
    adjusted_annots <- meth_annotated %>%
        group_by(chr, region) %>%
        left_join(first_regions) %>%
        mutate(offset_pos = ifelse(region == first_region,
                                   # no offset for first region
                                   pos,
                                   # offset for subsequent regions
                                   pos + cumsum_pos
                                   )) %>%
        ungroup() %>%
        arrange(chr, region, offset_pos) 
    adjusted_annots[40:50,]
    head(adjusted_annots)

    (nsig <- sum(adjusted_annots$fdr < fdr_thresh))

    meth_annotated <- adjusted_annots %>%
        select(chr, pos=offset_pos, stat, fdr, diff)


    # create granges for dmracte
    annotated <- GRanges(as.character(meth_annotated$chr),
                         IRanges(meth_annotated$pos, meth_annotated$pos),
                         stat = meth_annotated$stat,
                         diff = meth_annotated$diff,
                         ind.fdr = meth_annotated$fdr,
                         is.sig = meth_annotated$fdr < fdr_thresh
    )
    names(annotated) <- rownames(meth_annotated)
    annotated <- sort(annotated)
    head(annotated)

    head(meth_annotated)

    cpg_annotated <- new("CpGannotated", ranges=annotated)

    dmr_results <- tryCatch(
             {
                dmrcate(cpg_annotated, 
                                       lambda=1000, 
                                       C=2,
                                       min.cpgs=2
                                    )
             },
             error = function(cond){
                 message('error occurred')
                 message(cond$message)
                 return(NULL)
             },
             warning = function(cond){
                 message('warning occurred')
                 message(cond$message)
                 return(NULL)
             }
    )
    dmr_results
    #str(dmr_results)

    if (all(is.na(dmr_results))) return(NULL)

    # get the results as a dataframe
    dmr_df <- data.frame(coord = dmr_results@coord,
                         no.cpgs = dmr_results@no.cpgs,
                         min_smoothed_fdr = dmr_results@min_smoothed_fdr,
                         Stouffer = dmr_results@Stouffer,
                         HMFDR = dmr_results@HMFDR,
                         Fisher = dmr_results@Fisher,
                         maxdiff = dmr_results@maxdiff,
                         meandiff = dmr_results@meandiff
    )
    head(dmr_df)
    nrow(dmr_df) # number of DMRs found

    dmr_res <- list(dmr_res=dmr_df, annots=adjusted_annots) 
    dmr_res
}

add_label <- function(df){
    df <- df %>%
        mutate(label = NA,
               label =case_when(is_dm & is_sig ~ 'tp',
                              !is_dm & is_sig ~ 'fp',
                              !is_dm & !is_sig ~ 'tn',
                              is_dm & !is_sig ~ 'fn'
               )
    ) 
    df
}

summarize_metrics <- function(formatted_dm_results, save_runname){
    # plot the tp,tn,fp,fn metrics
    names(formatted_dm_results)

    # bring all results to the gene-level

    # cpgs
    cpgs_dm <- formatted_dm_results[['cpgs']] %>%
        select(feature, is_dm, is_sig) %>%
        separate(feature, c('chr', 'gene', 'cpg'), sep='_') %>%
        group_by(chr, gene) %>%
        summarize(is_dm = max(is_dm),
                  is_sig = max(is_sig)
        ) %>%
        unite('feature', chr, gene) %>%
        mutate(summary_type = 'cpgs')
    cpgs <- add_label(cpgs_dm) %>%
        select(feature, label, summary_type)
    head(cpgs)

    # use this as a full list of features
    cpgs_dm <- cpgs_dm %>%
        select(feature, is_dm) %>%
        unique()
    head(cpgs_dm)
    dim(cpgs_dm)

    # avgs already at gene-level
    avgs <- formatted_dm_results[['avgs']] %>%
        select(feature, is_dm, is_sig) %>%
        right_join(cpgs_dm) %>%
        replace_na(list(is_sig=0)) %>%
        mutate(summary_type = 'avgs')
    avgs <- add_label(avgs) %>%
        select(feature, label, summary_type)
    head(avgs)

    # make sure not duplicate pcs
    rpcs <- formatted_dm_results[['rpcs']] %>%
        select(feature, is_dm, is_sig) %>%
        separate(feature, c('chr', 'gene', 'pc'), sep='_') %>%
        group_by(chr, gene) %>%
        summarize(is_dm = max(is_dm),
                  is_sig = max(is_sig)
        ) %>%
        unite('feature', chr, gene) %>%
        right_join(cpgs_dm) %>%
        replace_na(list(is_sig=0)) %>%
        mutate(summary_type = 'rpcs')
    rpcs <- add_label(rpcs) %>%
        select(feature, label, summary_type)
    head(rpcs)


    # dmrs
    dmrs <- formatted_dm_results[['dmrs']] %>%
        select(feature, is_dm, is_sig) %>%
        separate(feature, c('chr', 'gene', 'cpg'), sep='_') %>%
        group_by(chr, gene) %>%
        summarize(is_dm = max(is_dm),
                  is_sig = max(is_sig)
        ) %>%
        unite('feature', chr, gene) %>%
        right_join(cpgs_dm) %>%
        replace_na(list(is_sig=0)) %>%
        mutate(summary_type = 'dmrs')
    dmrs <- add_label(dmrs) %>%
        select(feature, label, summary_type)
    head(dmrs)

    # combine the results together
    combined_res <- avgs %>%
        rbind(rpcs) %>%
        rbind(cpgs) %>%
        rbind(dmrs)
    head(combined_res)

    summarized_res <- combined_res %>%
        count(summary_type, label) %>%
        spread(label, n, fill=0) %>%
        mutate(precision = tp / (tp + fp),
               recall = tp / (tp + fn),
               specificity = tn / (tn + fp),
               sensitivity = tp / (tp + fn),
               accuracy = (tp + tn) / (tp + tn + fp + fn),

        )
    head(summarized_res)

    # plot the results
    #p <- ggplot(combined_res) +
        #geom_bar(aes(x=label,fill=summary_type),
                 #position=position_dodge()
        #)

    #(savefile <- paste0(savedir,
                        #"metrics_barplot",
                        #save_runname,
                        #".png"
        #))
    #ggsave(p, file=savefile)

    summarized_res
}

check_results <- function(dm_results , formatted_meth, fdr_thresh, save_runname){
    # check the DM results
    names(dm_results)
    names(formatted_meth)

    st <- 'rpcs'
    format_limma_results <- function(st){
        # get formatted input data
        df <- formatted_meth[[st]]
        meth <- df[['meth']]
        params <- df[['params']]
        head(params)

        # get output data
        dm_res <- dm_results[[st]]
        head(dm_res)

        # connect parameters to dm output 
        # to label tp, fp, tn, fn
        params$feature <- rownames(params)
        head(params)

        formatted_res <- dm_res %>%
            rownames_to_column('feature') %>%
            left_join(params, by='feature') %>%
            mutate(is_dm = as.numeric(percent_meth_difference > 0),
                   is_sig = as.numeric(bh_pval < fdr_thresh)
            ) %>%
            mutate(summary_type = st)
        formatted_res <- add_label(formatted_res)
        head(formatted_res)
        formatted_res
    }
    avgs <- format_limma_results('avgs')
    rpcs <- format_limma_results('rpcs')
    cpgs <- format_limma_results('cpgs')
    head(cpgs)

    # format the dmr results
    names(dm_results)
    names(formatted_meth)
    format_dmr_results <- function(){
        # ge the dmr results
        res <- dm_results[['dmrs']]
        names(res)


        # load true dmr ranges
        params <- formatted_meth[['cpgs']]$params
        params$feature <- rownames(params)
        head(params)

        # format dmr res to conne t annots
        dmr_res <- res[['dmr_res']] %>%
            separate(coord, c('chr', 'pos_range'), sep=':', remove=FALSE) %>%
            separate(pos_range, c('start', 'end'), sep='-')
        head(dmr_res)

        # get annoations for offset positions
        annots <- res[['annots']] %>%
            mutate(start = offset_pos,
                   end = offset_pos,
                   feature = paste(chr, region, pos, sep='_')
            ) %>%
            select(chr, start, end, feature)
        head(annots)

        # create genomicranges and find overlaps
        dmr_res_gr <- makeGRangesFromDataFrame(dmr_res, keep.extra.columns = TRUE)
        annots_gr <- makeGRangesFromDataFrame(annots, keep.extra.columns = TRUE)

        overlaps <- findOverlaps(subject=dmr_res_gr, query=annots_gr) %>%
            as.data.frame()

        dmr_res_annots <- cbind(as.data.frame(dmr_res_gr)[overlaps$subjectHits,],
                                annots=as.data.frame(annots_gr)[overlaps$queryHits,]
        ) %>%
            rename(feature = annots.feature) %>%
            select(!starts_with('annots'))
        head(dmr_res_annots)

        # connect the parameters
        head(params)
        added_params <- dmr_res_annots %>%
            left_join(params)
        head(added_params)

        # add labels
        labeled_dmrs <- added_params %>%
            mutate(is_dm = as.numeric(percent_meth_difference > 0),
                   is_sig = as.numeric(HMFDR < fdr_thresh)
            ) %>%
            mutate(summary_type = 'dmrs')
        labeled_dmrs <- add_label(labeled_dmrs)
        head(labeled_dmrs)
        labeled_dmrs
    }
    dmrs <- format_dmr_results()
    head(dmrs)

    formatted_dm_results <- list(avgs=avgs, rpcs=rpcs, cpgs=cpgs, dmrs=dmrs)
    summarized_metrics <- summarize_metrics(formatted_dm_results)

    #plot_metrics(formatted_dm_results, save_runname)

    summarized_dm <- list(formatted_dm_results=formatted_dm_results, 
                          summarized_dm_results=summarized_metrics)
    summarized_dm
}

plot_results <- function(){
    # decide how many genes
    total_genes <- 16000
    # how many genes have true DM?
    true_dm_proportion <- 0.05

    # need to fix the number of samples
    # to do DM on the matrix
    num_samples <- 50

    # repeat this run so that we can generate error bars on our final metrics
    repeat_times <- 100
    
    repeat_i <- 1
    run_analysis <- function(repeat_i){
        (save_runname <- paste0("_totalGenes", total_genes,
                               "_trueDmProportion", true_dm_proportion,
                               "_numSamples", num_samples,
                               "repeat", repeat_i
        ))
        message(paste("Loading results", repeat_i, "out of", repeat_times, Sys.time()))

        # skip if last output file exists
        (subsavedir <- paste0(savedir, "summarized_dm_res/"))
        (savefile <- paste0(subsavedir,
                           "summarized_dm_res",
                           save_runname,
                           ".rds"
        ))

        if (!file.exists(savefile)){
            message(paste("output file does not exist"))
            return(NA)
        }

        # read the results in
        res <- readRDS(savefile)
        names(res)

        # look at the summarized results
        summarized_res <- res[['summarized_dm_results']]
        head(summarized_res)
    }

}

main <- function(runnum){
   # meth simulation parameters
    # num_sites = number of CpG sites to simulate in the region
    # num_samples = number of samples to simulate data for
    # percent_meth_difference = how different should controls and cases be 
    # - at DMR sites
    # dmr_length = length of differentially methylated regions 
    # - this is not the same as # of differentially methylated CpGs
    # percent_sites_dm = percentages of CpGs that should be differentially
    # methylated 

    # ranges for parameters tested:
    # every combination of parameters was run
    #  num_sites_range <- c(20,50)
    #  num_samples_range <- c(50,500,5000)
    #  percent_meth_difference_range <- seq(0.1, 0.9, 0.1)
    #  percent_sites_dm_range <- seq(0,1,0.25)

    # parameters that do not change
    # N = 1000 for everything (number of regions simulated)
    # dmr_length = num_sites / 2

    # decide how many genes
    total_genes <- 16000
    #total_genes <- 100
    # how many genes have true DM?
    true_dm_proportion <- 0.05

    # need to fix the number of samples
    # to do DM on the matrix
    num_samples <- 50

    # repeat this run so that we can generate error bars on our final metrics
    repeat_times <- 100

    repeat_i <- 1
    run_analysis <- function(repeat_i){
        (save_runname <- paste0("_totalGenes", total_genes,
                               "_trueDmProportion", true_dm_proportion,
                               "_numSamples", num_samples,
                               "repeat", repeat_i
        ))
        message(paste("Running repeat analysis", repeat_i, "out of", repeat_times, Sys.time()))
        starttime <- Sys.time()

        # skip if last output file exists
        (subsavedir <- paste0(savedir, "summarized_dm_res/"))
        dir.create(subsavedir, showWarnings = FALSE)
        (savefile <- paste0(subsavedir,
                           "summarized_dm_res",
                           save_runname,
                           ".rds"
        ))
        if (file.exists(savefile)){
            message("Output file exists. Skipping run")
            return(1)
        }


        # create the data frame of mixed simulated data
        meth <- create_methylation_data(total_genes, true_dm_proportion, num_samples)
        head(meth$avgs)

        (subsavedir <- paste0(savedir, "meth_notNormalized/"))
        dir.create(subsavedir, showWarnings = FALSE)
        (savefile <- paste0(subsavedir,
                           "meth_notNormalized",
                           save_runname,
                           ".rds"
        ))
        saveRDS(meth, savefile)

        # normalize the methylation for DM analysis
        normalized_meth <- normalize_meth(meth)

        # format methylation for DM analysis
        formatted_meth <- format_meth(normalized_meth)

        (subsavedir <- paste0(savedir, "formatted_data/"))
        dir.create(subsavedir, showWarnings = FALSE)
        (savefile <- paste0(subsavedir,
                           "formatted_data",
                           save_runname,
                           ".rds"
        ))
        saveRDS(formatted_meth, savefile)
        #formatted_meth <- readRDS(savefile)

        # run DM analysis
        print('avgs')
        avgs_res <- run_dm_analysis(formatted_meth[['avgs']]$meth)
        print('rpcs')
        rpcs_res <- run_dm_analysis(formatted_meth[['rpcs']]$meth)
        print('cpgs')
        cpgs_res <- run_dm_analysis(formatted_meth[['cpgs']]$meth)

        # run dmr analysis - TO DO
        fdr_thresh <- 0.05
        dmr_res <- run_dmr_analysis(cpgs_res, fdr_thresh)
        names(dmr_res)
        head(dmr_res$dmr_res)
        dim(dmr_res$dmr_res)

        (subsavedir <- paste0(savedir, "dm_res/"))
        dir.create(subsavedir, showWarnings = FALSE)
        (savefile <- paste0(subsavedir,
                           "dm_res",
                           save_runname,
                           ".rds"
        ))
        dm_results <- list(avgs=avgs_res, rpcs=rpcs_res, cpgs=cpgs_res, dmrs=dmr_res)
        saveRDS(dm_results, savefile)

        # check the dm results
        summarized_dm <- check_results(dm_results, formatted_meth, fdr_thresh, save_runname)

        (subsavedir <- paste0(savedir, "summarized_dm_res/"))
        dir.create(subsavedir, showWarnings = FALSE)
        (savefile <- paste0(subsavedir,
                           "summarized_dm_res",
                           save_runname,
                           ".rds"
        ))
        saveRDS(summarized_dm, savefile)

        endtime <- Sys.time()
        duration <- difftime(endtime, starttime, units='secs')[[1]] / 60
        message(paste("Run", repeat_i, "out of", repeat_times, "completed in", round(duration, 2), "minutes"))
        1
    }

    set.seed(runnum)
    #res <- lapply(1:repeat_times, run_analysis)
    run_analysis(runnum)
}



# get command line arguments and run main
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
# run number
runnum = args[1]

main(runnum)
