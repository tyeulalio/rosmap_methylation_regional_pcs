library(RnBeads)
library(ppcor)
library(lmtest)

# for plotting
library(ggplot2)
library(RColorBrewer)
library(UpSetR)
library(ggridges)

library(tidyverse)


datadir <- paste0("/path/to/data/")
savedir <- paste0("/path/to/output/")


load_pheno <- function(svddir, cell_type){
    # load the phenotype data   
    datafile <- paste0(svddir, "astro_formatted_clinical.rds")
    pheno <- readRDS(datafile)
    head(pheno)

    # load global meth PCs also
    datafile <- paste0(svddir, "cleaned_", cell_type, "_global_pcs.rds")
    global_pcs <- readRDS(datafile) %>%
        rownames_to_column('individualID')
    head(global_pcs)

    # set up variables 
    sub_pheno <- pheno %>%
        select(individualID, braaksc, ceradsc, cogdx, dcfdx_lv, diagnosis, apoe4,
               age_death, pmi, sex, batch=meth.batch, neuron:endo, Study
        ) %>%
        mutate(Study = as.numeric(as.factor(Study))) %>%
        left_join(global_pcs)
    head(sub_pheno)

    # select the numeric and factor variables
    numeric_vars = c('braaksc', 'ceradsc', 'cogdx', 'dcfdx_lv', 'apoe4', 'age_death', 'pmi', 'Study', 'neuron', 'oligo_opc', 'astro', 'endo')
    factor_vars = c('diagnosis', 'sex', 'batch')
    factor_nums = paste0("cat_", numeric_vars)
    
    # transform variables here
    sub_pheno[,numeric_vars] = lapply(sub_pheno[,numeric_vars], as.character)
    sub_pheno[,numeric_vars] = lapply(sub_pheno[,numeric_vars], as.numeric)
    sub_pheno[,factor_vars] = lapply(sub_pheno[,factor_vars], as.factor)

    # test categorical numeric vars also
    sub_pheno[,factor_nums] = lapply(sub_pheno[,numeric_vars], as.factor)

    # remove any if age_Death and sex are missing
    sub_pheno <- sub_pheno %>%
        filter(!is.na(age_death), !is.na(sex))

    sub_pheno
}


load_meth <- function(region_type, summary_type, cell_type){
    # load the methylation data

    # get the data file depending on summary type
    datafile <- NA
    if (summary_type == 'cpgs'){
        datafile <- paste0(datadir, cleaned_str, cell_type, "_lifted_hg38_mvals.rds")
    } else{
        datadir
        list.files(datadir)

        datafile <- paste0(datadir, cleaned_str, region_type, "_", cell_type, "_", summary_type, ".rds")
    }

    # load the datafile
    datafile
    meth <- readRDS(datafile)
    head(meth)
    dim(meth)


    # update datafile so they all match depending on the summary type
    if (summary_type == 'avgs'){
        meth <- meth %>%
            select(-duration_secs)
    }
    if (summary_type == 'pcs'){
        # remove columns that are not pcs
        head(meth)
        meth <- meth %>%
            select(-percent_variance_explained)
    }

    # load cpg map also
    datafile <- paste0(datadir, region_type, "_", cell_type, "_gene_cpg_map.rds")
    cpg_map <- readRDS(datafile)
    head(cpg_map)
    
    list(meth=meth, cpg_map=cpg_map)
}

count_regions <- function(runs){
    # count the cpgs, and pcs in each region
    runs <- expand.grid(region_type=region_types,
                        cell_type=cell_types
    )
    head(runs)
    dim(runs)

    row <- runs[1,]
    load_data <- function(row){
        region_type <- row[['region_type']] %>% as.character()
        cell_type <- row[['cell_type']] %>% as.character()

        print(paste("Processing", region_type, cleaned_str, cell_type))

        summary_type = 'pcs'

        dat <- load_meth(region_type, summary_type, cell_type)
        names(dat)


        # check the genes and number of pcs per gene
        meth <- dat$meth
        cpg_map <- dat$cpg_map

        # count cpgs first
        head(meth)
        gene_pcs <- data.frame(gene_pc=rownames(meth)) %>%
            separate(gene_pc, c('gene_id', 'pc'), sep='-') 

        # get cpg-level stats from cpg_map, including
        # total genes, total CpGs, CpGs per gene
        head(cpg_map)
        overall_stats <- cpg_map %>%
            summarize(total_genes = n_distinct(gene_id),
                      total_cpgs = n_distinct(cpg_id),
            ) 
        head(overall_stats)

        # cpgs per gene
        cpg_stats <- cpg_map %>%
            group_by(gene_id) %>%
            summarize(cpg_count = n_distinct(cpg_id)) %>%
            ungroup() %>%
            summarize(min_cpg = min(cpg_count),
                      mean_cpg = mean(cpg_count),
                      median_cpg = median(cpg_count),
                      max_cpg = max(cpg_count)
            ) 
        head(cpg_stats)

        # get pc-level stats
        # number of pcs per gene, number of PCs total
        head(gene_pcs)
        length(unique(gene_pcs$gene_id))
        pc_stats <- gene_pcs %>%
            group_by(gene_id) %>%
            summarize(pc_count = n_distinct(pc)) %>%
            ungroup() %>%
            summarize(min_pc = min(pc_count),
                      mean_pc = mean(pc_count),
                      median_pc = median(pc_count),
                      max_pc = max(pc_count)
            ) 
        head(pc_stats)


        # load variance explained
        list.files(datadir)
        datafile <- paste0(datadir, cleaned_str, region_type, "_", cell_type, "_pcs.rds")
        var_explained <- readRDS(datafile)[,1,drop=FALSE]
        head(var_explained)

        var_stats <- var_explained %>%
            summarize(min_var_expl = min(percent_variance_explained),
                      mean_var_expl = mean(percent_variance_explained),
                      median_var_expl = median(percent_variance_explained),
                      max_var_expl = max(percent_variance_explained)
            ) 
        var_stats

        combined_stats = cbind(overall_stats, cpg_stats) %>%
            cbind(pc_stats) %>%
            cbind(var_stats) %>%
            mutate(region_type = region_type,
                   cell_type = cell_type
            ) %>%
            select(region_type, cell_type, everything())
        head(combined_stats)

        combined_stats
    }

    celltype_dat <- apply(runs, 1, load_data) 
    counts_df <- do.call(rbind, celltype_dat)
    head(counts_df)

    (savefile <- paste0(savedir, "feature_per_gene_counts.csv"))
    write.csv(counts_df, savefile)
}


correlate_gene_complexity <- function(runs){
    # correlate number of rPCs with the gene complexisty
    runs <- expand.grid(region_type=region_types,
                        cell_type=cell_types
    )
    head(runs)
    dim(runs)

    row <- runs[1,]
    load_data <- function(row){
        region_type <- row[['region_type']] %>% as.character()
        cell_type <- row[['cell_type']] %>% as.character()

        print(paste("Processing", region_type, cleaned_str, cell_type))

        summary_type = 'pcs'

        dat <- load_meth(region_type, summary_type, cell_type)
        names(dat)


        # check the genes and number of pcs per gene
        meth <- dat$meth
        cpg_map <- dat$cpg_map

        # how many cpgs per promoter?
        head(cpg_map)
        cpg_counts <- cpg_map %>%
            count(gene_id)
        head(cpg_counts)

        # look at pcs per gene
        gene_pcs <- data.frame(gene_pc=rownames(meth)) %>%
            separate(gene_pc, c('gene_id', 'pc'), sep='-') 
        head(gene_pcs)
        pc_counts <- gene_pcs %>%
            count(gene_id) 
        head(pc_counts)

        # relationship between CpGs and PCs per gene
        head(gene_pcs)
        head(cpg_counts)

        pc_cpgs <- gene_pcs %>%
            count(gene_id, name='num_pcs') %>%
            inner_join(cpg_counts) %>%
            rename(num_cpgs=n) %>%
            mutate(cell_type=cell_type,
                   region_type=region_type
            )
        head(pc_cpgs)

        pc_cpgs
    }

    celltype_dat <- apply(runs, 1, load_data) 
    pc_cpgs <- do.call(rbind, celltype_dat)
    head(pc_cpgs)

    (savefile <- paste0(savedir, "gene_pc_cpgs_joined.csv"))
    write_csv(pc_cpgs, savefile)

    plotdir <- paste0(savedir, "plots/")
    dir.create(plotdir)
    plotdir

    head(pc_cpgs)
    tail(pc_cpgs)

    unique(pc_cpgs$cell_type)
    p <- pc_cpgs %>%
        left_join(cell_type_map) %>%
        left_join(region_type_map) %>%
        #filter(region_type == 'full_gene') %>%
        ggplot(aes(x=num_cpgs, y=num_pcs)) +
        geom_point(aes(color=formatted_ct), size=3) +
        geom_smooth(aes(color=formatted_ct), method='loess', se=FALSE) +
        facet_wrap(. ~ formatted_rt, nrow=2, scales='free') +
        scale_color_manual(values=cell_type_colors) +
        theme_bw() +
        theme(text=element_text(size=20),
              ) +
        common_theme +
        guides(color=guide_legend(title=paste("Cell type"))) +
        xlab(paste("CpGs per gene region")) +
        ylab(paste("rPCs per gene region")) 
    #(savefile <- paste0(plotdir, "pc_cpg_scatter.png"))
    (savefile <- paste0(plotdir, "pc_cpg_scatter_all_regions.png"))
    ggsave(p, file=savefile, width=10, height=10, unit='in')


    # get correlations for each cell type separately
    ct <- 'bulk'
    rt <- 'promoters'
    check_cell_type <- function(row){
        ct <- row[['cell_type']]
        rt <- row[['region_type']]
        print(paste("Processing", ct, rt))
        # what is correlation between CpGs and PCs?
        #cor.test(pc_cpgs$num_pcs, pc_cpgs$num_cpgs, method='spearman')

        # use an interaction model to summarize all cell types
        #lm_res <- lm(num_pcs ~ cell_type * num_cpgs + 0, data=pc_cpgs)
        #summary(lm_res)

        ct_pc_cpgs <- pc_cpgs %>%
            filter(cell_type == ct,
                   region_type == rt
            )
        head(ct_pc_cpgs)
        unique(ct_pc_cpgs$num_pcs)

        cpg_corr <- cor.test(ct_pc_cpgs$num_pcs, ct_pc_cpgs$num_cpgs, method='pearson')

        # compute coeffieint of variation for each gene to plot against num PCs
        # load mvals
        datafile <- paste0(datadir, cleaned_str, ct, "_lifted_hg38_mvals.rds")
        mvals <- readRDS(datafile)
        head(mvals)

        # for each gene, compute the coefficient of variantion
        summary_type = 'pcs'
        dat <- load_meth(rt, summary_type, ct)
        names(dat)

        # check the genes and number of pcs per gene
        cpg_map <- dat$cpg_map
        head(cpg_map)

        #head(mvals)
        #row_vars <- apply(mvals,1,sd)
        #head(row_vars)

        ct_cpg_map <- cpg_map 

        head(ct_cpg_map)
        gene_id <- ct_cpg_map[1,]$gene_id
        compute_variation <- function(gene_id){
            # get cpgs for this gene
            gene_cpgs <- ct_cpg_map[which(ct_cpg_map$gene_id == gene_id),]$cpg_id
            head(gene_cpgs)

            # get mvals for cpgs
            gene_vars <- mvals[gene_cpgs,]
            head(gene_vars)

            #vars <- mean(gene_vars)
            vars = mad(gene_vars)
            
            tibble(gene_id=gene_id, var=vars)
        }

        genes <- unique(ct_cpg_map$gene_id) 

        # takes a while to run
        (savefile <- paste0(savedir, ct, "_", rt, "_mads.rds"))
        if (file.exists(savefile)){
            print(paste("LOADING EXISTING MADS FILE"))
            mads_df <- readRDS(savefile)
        } else{
            Sys.time()
            mads <- lapply(genes, compute_variation)
            Sys.time()

            mads_df <- do.call(rbind, mads)
            head(mads_df)
            saveRDS(mads_df, savefile)
        }

        # find relationshp between variation and PCs
        pc_mads <- mads_df %>%
            left_join(ct_pc_cpgs)
        head(pc_mads)


        pc_mads_corr <- cor.test(pc_mads$var, pc_mads$num_pcs, method='pearson')
        #cor.test(pc_mads$var, pc_mads$num_cpgs, method='pearson')

        


        # get the length of genes and correlate with number of PCs
        #print("checking full gene length**")
        #region_type <- 'full_gene'
        datafile <- paste0(datadir, "annotated_regions/",
                           rt, "_annotations.rds"
        )
        annots <- readRDS(datafile)
        head(annots)

        if (rt == 'full_gene'){
            unique(annots$type)
            sub_annots <- annots %>%
                filter(type %in% c('hg38_genes_1to5kb', 'hg38_genes_3UTRs')) %>%
                mutate(type = str_remove(type, "hg38_genes_")) %>%
                select(seqnames, start, end, type, gene_id) %>%
                filter(end > start)
            head(sub_annots)
        } else{
            sub_annots <- annots %>%
                mutate(type = str_remove(type, "hg38_genes_")) %>%
                select(seqnames, start, end, type, gene_id) %>%
                filter(end > start)
        }
        head(sub_annots)

        wide_start <- sub_annots %>%
            group_by(gene_id) %>%
            summarize(start = min(start)) 
        head(wide_start)

        wide_end <- sub_annots %>%
            group_by(gene_id) %>%
            summarize(end = max(end)) 
        head(wide_start)

        head(pc_cpgs)
        gene_boundaries <- inner_join(wide_start, wide_end) %>%
            mutate(length = end-start) %>%
            left_join(ct_pc_cpgs) %>%
            filter(!is.na(num_pcs))
        head(gene_boundaries)

        summary(gene_boundaries$length)

        length_corr <- cor.test(gene_boundaries$length, gene_boundaries$num_pcs, method='pearson')

        data.frame(mads_corr=pc_mads_corr$estimate,
                   mads_pval=pc_mads_corr$p.value,
                   length_corr=length_corr$estimate,
                   length_pval=length_corr$p.value,
                   cpg_corr=cpg_corr$estimate,
                   cpg_pval=cpg_corr$p.value,
                   cell_type=ct,
                   region_type=rt
        )
    }

    all_res <- apply(runs, 1, check_cell_type)
    res_df <- do.call(rbind, all_res)
    head(res_df)

    (savefile <- paste0(savedir, "gene_complexity_correlations.csv"))
    write_csv(res_df, savefile)


    # plot the number of features
    head(cpg_map)
    head(mvals)
    head(gene_pcs)

    # how many cpgs?
    n_cpgs <- length(unique(cpg_map$cpg_id))

    # how many pcs?
    n_pcs <- dim(gene_pcs)[1]

    # how many avgs?
    n_avgs <- length(unique(gene_pcs$gene))

    counts_df <- data.frame(summary_type=c('CpGs', 'PCs', 'Avgs'), count=c(n_cpgs, n_pcs, n_avgs))
    counts_df    


    p <- ggplot(counts_df) +
        geom_bar(aes(x=summary_type, y=count, fill=summary_type),
                 stat='identity',
                 show.legend=FALSE
        ) +
        geom_label(aes(x=summary_type, y=count, label=count),
                   size=8) +
        scale_fill_manual(values=summary_type_colors) +
        theme_bw() +
        theme(text=element_text(size=20)) +
        xlab(paste("Summary type")) +
        ylab(paste("Total features"))
    savefile <- paste0(plotdir, "feature_count.png")
    ggsave(p, file=savefile)


    # check percent variance percent_variance_explained
    head(runs)
    check_pct_var <- function(row){
        region_type <- row[['region_type']] %>% as.character()
        summary_type <- row[['summary_type']] %>% as.character()
        cell_type <- row[['cell_type']] %>% as.character()

        print(paste("Processing", region_type, summary_type, cell_type))

        list.files(datadir)
        filename <- paste(region_type, cell_type, "pcs.rds", sep='_')
        (datafile <- paste0(datadir, cleaned_str, filename))
        pc_info <- readRDS(datafile)[1]
        head(pc_info)

        pc_info %>%
            rename(x=percent_variance_explained) %>%
            summarise(min=min(x),
                      mean=mean(x),
                      median=median(x),
                      max=max(x),
                      cell_type=cell_type
            )
    }
    res <- apply(runs,1,check_pct_var)
    res_df <- do.call(rbind, res)
    res_df

    savefile <- paste0(savedir, "percent_variance_explained.csv")
    write_csv(res_df, savefile)

}

get_global_pcs <- function(cell_type, region_type){
    # run pca
    list.files(svddir)
    datafile <- paste0(svddir, cleaned_str, "cleaned_", cell_type, "_global_pcs.rds")
    datafile
    global_pcs <- readRDS(datafile)
    head(global_pcs)

    global_pcs[-1,]
}

correlate_global_local_pcs <- function(meth, global_pcs, cell_type, region_type, summary_type){
    # correlate local pcs with global pcs
    head(meth)
    head(global_pcs)

    t_global <- t(global_pcs)
    rownames(t_global) <- paste0("global_", rownames(t_global))
    head(t_global)

    # make sure samples are in order
    ordered_samples <- intersect(colnames(t_global), colnames(meth))
    head(ordered_samples)
    length(ordered_samples)

    ordered_local <- meth[,ordered_samples]
    ordered_global <- t_global[,ordered_samples]

    if (!identical(colnames(ordered_local), colnames(ordered_global))){
        print("ORDER DOESN'T MATCH")
        return()
    }

    # get corrs from the ordered mattices
    corrs <- cor(t(ordered_local), t(ordered_global), method='pearson')
    head(corrs)

    format_corrs <- function(corrs){
        formatted_corrs <- corrs %>%
            as.data.frame() %>%
            rownames_to_column('gene_pc') %>%
            separate(gene_pc, c('gene_id', 'pc'), sep='-') %>%
            gather('global', 'corr', -gene_id, -pc) %>%
            mutate(corr = abs(corr))
        head(formatted_corrs)

        # get the maximum corr for each gene
        max_corrs <- formatted_corrs %>%
            group_by(gene_id) %>%
            summarise(max_corr = max(corr))
        head(max_corrs)

        list(formatted_corrs=formatted_corrs,
             max_corrs=max_corrs
        )
    }

    # get regular correlations
    res <- format_corrs(corrs)
    formatted_corrs <- res$formatted_corrs
    max_corrs <- res$max_corrs

    plotdir <- paste0(savedir, "plots/")
    dir.create(plotdir)

    plot_res <- function(formatted_corrs, max_corrs, random){
        head(max_corrs)
        (median <- mean(max_corrs$max_corr))


        # check max corr split by pc number
        head(formatted_corrs)

        pc_max_corrs <- formatted_corrs %>%
            group_by(gene_id, pc) %>%
            summarise(max_corr = max(corr))
        head(pc_max_corrs)

        filter_pcs <- pc_max_corrs %>%
            ungroup() %>%
            count(pc) %>%
            filter(n < 2)
        filter_pcs

        pc_max_corrs <- pc_max_corrs %>%
            filter(!pc %in% filter_pcs$pc)


        pc_max_corrs$pc <- factor(pc_max_corrs$pc,
                               levels=paste0("PC", 1:length(unique(pc_max_corrs$pc))
        ))
        
        # plot two distributions together
        head(max_corrs)
        head(formatted_corrs)

        maxc <- max_corrs %>%
            rename(corr=max_corr) %>%
            mutate(corr_type = 'Max')
        head(maxc)

        formc <- formatted_corrs %>%
            select(-pc, -global) %>%
            mutate(corr_type = 'Overall')
        head(formc)


        joined_corrs <- rbind(maxc, formc)
        head(joined_corrs)

        # get median for each
        meds <- joined_corrs %>%
            group_by(corr_type) %>%
            summarize(med=median(corr))
        meds

        p <- ggplot(joined_corrs) +
            geom_density(aes(x=corr, fill=corr_type),
                         alpha=0.6) +
            geom_vline(xintercept=meds[[1,'med']], color='red') +
            geom_vline(xintercept=meds[[2,'med']], color='blue') +
            geom_label(data=meds, aes(x=med, y=9, color=corr_type, label=round(med,2)),
                       show.legend=FALSE) +
            coord_cartesian(xlim=c(0,1)) +
            theme_bw() +
            theme(text=element_text(size=20)) +
            ggtitle(paste("Max local,global PC correlation", region_type, cell_type)) +
            guides(fill=guide_legend(title=paste("Correlation"))) +
            xlab(paste("Pearson correlation")) +
            ylab(paste("Density"))
        (savefile <- paste0(plotdir, summary_type, cell_type, "_overlapped_correlation_density.png"))
        ggsave(p, file=savefile)
        
        1
    }

    plot_res(formatted_corrs, max_corrs, random=FALSE)
    1
}

check_global_pcs <- function(dat, cell_type, region_type, summary_type){
    # check correlation of local pcs with global pcs
    names(dat)
    meth <- dat$meth
    cpg_map <- dat$cpg_map

    # get global pcs
    global_pcs <- get_global_pcs(cell_type, region_type)
    head(global_pcs)

    # correlate pcs with global pcs
    correlate_global_local_pcs(meth, global_pcs, cell_type, region_type, summary_type)
}

check_pc_correlations <- function(meth, pcs_cpg){
    # check correlations between PCs across genes
    head(meth)[1:10,]
    head(pcs_cpg)

    # get a map from gene to chromosome
    chrom_map <- pcs_cpg %>%
        select(gene_id, pc, cpg_id) %>%
        separate(cpg_id, c('chrom', 'start', 'end'), sep='_') %>%
        select(-start, -end) %>%
        mutate(gene_pc = paste(gene_id, pc, sep='-')) %>%
        unique()
    head(chrom_map)

    global_pc_num <- 'PC0'
    check_global_pcs <- function(global_pc_num){

    # test PCs
    pc_num <- 'PC3'
    test_pc <- function(pc_num){
        #(savefile <- paste0(savedir, cleaned_str, pc_num, "_pc_partial_correlation_with", global_pc_num, ".rds"))
        #if (file.exists(savefile)){
            #print("Loading existing results.")
            #res_df <- readRDS(savefile)
            #return(res_df)
        #}

        pc_genes <- chrom_map %>%
            filter(pc == pc_num) 
        head(pc_genes)
        dim(pc_genes)
        table(pc_genes$chrom)

        # subset PCs
        set.seed(1174117174)
        sub_genes <- pc_genes[sample(1:nrow(pc_genes), 100, replace=FALSE),]
        pc_genes <- sub_genes

        table(pc_genes$chrom)

        comparisons <- expand.grid(gene1 = pc_genes$gene_pc, gene2 = pc_genes$gene_pc)
        head(comparisons)
        dim(comparisons)

        # only look at comparisons across chromosomes
        filtered_comparisons <- comparisons %>%
            # get chromosome for gene1
            mutate(gene_pc = gene1) %>%
            left_join(pc_genes[c('gene_pc', 'chrom')]) %>%
            select(-gene_pc) %>%
            rename(chrom1 = chrom) %>%
            # get chromosomes for gene2
            mutate(gene_pc = gene2) %>%
            left_join(pc_genes[c('gene_pc', 'chrom')]) %>%
            select(-gene_pc) %>%
            rename(chrom2 = chrom) %>%
            # remove matched chromosomes
            filter(chrom1 != chrom2) %>%
            mutate(rn = row_number())
        head(filtered_comparisons)
        dim(filtered_comparisons)
        rm(comparisons)

        # subset methylation to this level
        head(meth)
        sub_meth <- meth[pc_genes$gene_pc,]
        dim(sub_meth)
        head(sub_meth)

        # add phenotype data that we account for
        head(pheno)
        sub_pheno <- pheno %>%
            #select(individualID, age_death, pmi, sex, batch, Study, neuron:endo, PC1:PC27) %>%
            select(individualID, age_death, pmi, sex, batch, Study, neuron:endo) %>%
            mutate(sex = as.numeric(as.factor(sex))) %>%
            column_to_rownames('individualID') 
        head(sub_pheno)

        sub_pheno <- sub_pheno[!colnames(sub_pheno) %in% c(global_pc_num)]
        head(sub_pheno)

        ordered_pheno <- sub_pheno[colnames(sub_meth),]
        identical(rownames(ordered_pheno), colnames(sub_meth))
        head(ordered_pheno)

        row <- filtered_comparisons[1,]
        check_pair <- function(row){
            gene1 <- row[['gene1']] %>% as.character()
            gene2 <- row[['gene2']] %>% as.character()
            rn <- row[['rn']] %>% as.numeric()

            if (rn %% 1000 == 0) print(paste("Processed row", rn, "out of", nrow(filtered_comparisons), Sys.time()))

            gene_dat1 <- sub_meth[c(gene1, gene2),] %>%
                t() %>%
                as.data.frame() %>%
                rownames_to_column('individualID')
            dim(gene_dat1)
            head(gene_dat1)

            ordered_pheno['individualID'] <- rownames(ordered_pheno)

            head(ordered_pheno)
            dim(ordered_pheno)
            pair_dat <- gene_dat1 %>%
                #left_join(ordered_pheno, by='individualID') %>%
                select(-individualID) 
                #mutate(pmi = ifelse(is.na(pmi), mean(pmi, na.rm=TRUE), pmi)) 
            head(pair_dat)

            pair_mat <- data.matrix(pair_dat)
            head(pair_mat)

            partial_corr <- pcor(pair_mat, method='spearman')
            partial_corr

            i <- 1
            j <- 2
            par_corr <- partial_corr$estimate[i,j]
            pval <- partial_corr$p.value[i,j]

            data.frame(gene1=gene1, gene2=gene2, partial_corr=par_corr, pval=pval)
        }
        res <- apply(filtered_comparisons, 1, check_pair)
        res_df <- do.call(rbind, res) %>%
            mutate(pc = pc_num)
        head(res_df)

        (savefile <- paste0(savedir, cleaned_str, pc_num, "_pc_partial_correlation_with", global_pc_num, ".rds"))
        saveRDS(res_df, savefile)

        res_df

    }

    pcs <- paste0('PC', 1:3)
    res <- lapply(pcs, test_pc)
    res_df <- do.call(rbind, res)
    head(res_df)

    formatted_res <- res_df

    (savefile <- paste0(savedir, "new_", cleaned_str, "pc_partial_correlation_distributions_with", global_pc_num, ".png"))
    p <- ggplot(formatted_res) +
        #geom_density(aes(x=abs(corr), color=pc), linewidth=2) +
        geom_density(aes(x=abs(partial_corr), color=pc), linewidth=2) +
        theme_bw() +
        theme(text=element_text(size=18))
    ggsave(p, file=savefile)
    1
    }

    check_global_pcs("PC0")

    global_pc_nums <- paste0("PC", 1:27)
    lapply(global_pc_nums, check_global_pcs) 


    # load the data
    pc_num <- 'PC1'
    combine_results <- function(pc_num){

        global_pc_num <- 'PC1'
        load_global <- function(global_pc_num){
            (datafile <- paste0(savedir, pc_num, "_pc_partial_correlation_with", global_pc_num, ".rds"))
            part_corrs <- readRDS(datafile) %>%
                mutate(global_pc = global_pc_num)
            part_corrs
        }
        global_pcs <- paste0("PC", 0:27)
        global_res <- lapply(global_pcs, load_global) %>%
            do.call(rbind, .) %>%
            mutate()
        head(global_res)

        global_res
    }
    local_pcs <- paste0("PC", 1:3)
    res <- lapply(local_pcs, combine_results) %>%
        do.call(rbind, .)
    head(res)

    global_pcs <- paste0("PC", 0:5)
    (savefile <- paste0(savedir, "new_combined_pc_partial_correlations.png"))
    p <- res %>%
        filter(global_pc %in% global_pcs) %>%
        mutate(global_pc = factor(global_pc, levels=global_pcs)) %>%
        ggplot() +
        geom_density(aes(x=abs(partial_corr), color=pc), linewidth=1.5) +
        facet_wrap(. ~ global_pc, scales='free_y', nrow=2) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=50, hjust=1, vjust=1),
              strip.background=element_rect(fill='white')
        )
    ggsave(p, file=savefile, width=8.5, height=4)
}

main <- function(){
    region_types <- c('full_gene', 'promoters', 'gene_body', 'preTSS')
    summary_types <- c('pcs', 'avgs')

    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')

    runs <- expand.grid(region_type=region_types,
                        summary_type=summary_types,
                        cleaned_str=cleaned_strs,
                        cell_type=cell_types
    )
    head(runs, 10)
    runs
    dim(runs)


    # count features per gene region
    count_regions(runs)

    # correlate gene complexity with rPCs
    correlate_gene_complexity(runs)
    
    row <- runs[1,]
    process_row <- function(row){
        region_type <- row[['region_type']] %>% as.character()
        summary_type <- row[['summary_type']] %>% as.character()
        cell_type <- row[['cell_type']] %>% as.character()

        if (region_type == 'full_array' & summary_type != 'cpgs') return()

        print(paste("Processing", region_type, summary_type, cleaned_str, cell_type))

        # load phenotype data
        pheno <- load_pheno(cell_type)
        head(pheno)


        # load methylation data
        dat <- load_meth(region_type, summary_type, cell_type)
        meth <- dat$meth
        cpg_map <- dat$cpg_map
        head(meth)
        head(meth)[1:10]

        # checking global pcs
        check_global_pcs(dat, cell_type, region_type, summary_type)

        return(data.frame(region_type=region_type,
                          summary_type=summary_type,
                          cell_type=cell_type,
                          num_genes=nrow(meth)))




        gene_pcs <- data.frame(gene_pcs=rownames(meth)) %>%
            separate(gene_pcs, c('gene_id', 'pc'), sep='-')
        head(gene_pcs)

        head(cpg_map)

        pcs_cpg <- gene_pcs %>%
            full_join(cpg_map) %>%
            mutate(region_type = region_type,
                   cell_type = cell_type
            )
        head(pcs_cpg)

        check_pc_correlations(meth, pcs_cpg)


        1
    }

    res <- apply(runs, 1, process_row)
    res_df <- do.call(rbind, res)
    res_df
    
}

main()
