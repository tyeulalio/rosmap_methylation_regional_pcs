library(limma)
library(RnBeads.hg38)
library(RnBeads)
library(RNOmni)
library(isva)
library(biomaRt)

library(lme4)
library(broom.mixed)

# for plotting
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(UpSetR)
library(ComplexUpset)
library(ggupset)

library(tidyverse)






# use command line input for cell type
args = commandArgs(trailingOnly=TRUE)
celltype_idx <- args[1] %>%
    as.numeric()

# define the input and output directories
datadir <- "/path/to/summarized/methylation/"
savedir <- "/path/to/save/"
dir.create(savedir, showWarnings=FALSE)




# include global pcs as covariates or not?
include_global = TRUE
protect_global = TRUE
if (include_global){
    savedir <- paste0(savedir, "global_pcs/")
    dir.create(savedir)

    if (protect_global){
        savedir <- paste0(savedir, "protect_global_pcs/")
        dir.create(savedir)
    } else{
        savedir <- paste0(savedir, "remove_all_global_pcs/")
        dir.create(savedir)
    }

}

cleaned_mvals = TRUE
if (cleaned_mvals){
    cleaned_str = "cleaned_"
} else{
    cleaned_str = ""
}
clean_global = TRUE
if (clean_global){
    cleaned_str = "cleaned_global_"
}

test_bulk = FALSE
if (test_bulk){
    savedir <- paste0(savedir, "test_bulk/")
    dir.create(savedir)
}

datadir
savedir


load_pheno <- function(cell_type){
    # load the phenotype data   
    # can use any cell type since they're all the same 
    datafile <- paste0(svddir, "astro_formatted_clinical.rds")
    pheno <- readRDS(datafile)
    head(pheno)

    # set up variables 
    sub_pheno <- pheno %>%
        select(individualID, braaksc, ceradsc, cogdx, dcfdx_lv, diagnosis, apoe4, apoe_genotype, bin_apoe4,
               age_death, pmi, sex, batch=meth.batch, study=Study,
               neuron, oligo_opc, astro, endo, new_apoe, new_diag
        ) %>%
        mutate(new_diag = ifelse(new_diag == 'Control', NA, new_diag),
                new_diag = as.numeric(new_diag == 'AD')
        ) %>%
        mutate(new_apoe = as.numeric(new_apoe == 34)) 
    head(sub_pheno)


    # format numeric and factor variables
    numeric_vars = c('braaksc', 'ceradsc', 'cogdx', 'dcfdx_lv', 'age_death', 'pmi',
                     'neuron', 'oligo_opc', 'astro', 'endo', 'apoe4', 'bin_apoe4',
                     'new_apoe', 'new_diag'
    )
    factor_vars = c('diagnosis', 'sex', 'batch', 'apoe_genotype')
    factor_nums = paste0("cat_", numeric_vars)
    
    sub_pheno[,numeric_vars] = lapply(sub_pheno[,numeric_vars], as.character)
    sub_pheno[,numeric_vars] = lapply(sub_pheno[,numeric_vars], as.numeric)
    sub_pheno[,factor_vars] = lapply(sub_pheno[,factor_vars], as.factor)

    # test categorical numeric vars also
    sub_pheno[,factor_nums] = lapply(sub_pheno[,numeric_vars], as.factor)

    # remove any if age_Death and sex are missing
    sub_pheno <- sub_pheno %>%
        filter(!is.na(age_death), !is.na(sex))
    head(sub_pheno)

    # add categorical for braak
    sub_pheno <- sub_pheno %>%
        mutate(cat_braaksc = as.numeric(braaksc >= 5)) %>%
        mutate(apoe_genotype = as.factor(apoe_genotype))
    head(sub_pheno)

    sub_pheno
}


run_dmp <- function(meth, pheno, global_pcs, region_type, summary_type, cleaned_str, trait, cell_type, shuffle_pheno, shuffle_num){
    # run the dmp analysis

    # attach global pcs to ordered pheno
    head(global_pcs)
    pc_cols = colnames(global_pcs)
    pc_cols

    formatted_pcs <- global_pcs %>%
        as.data.frame() %>%
        rownames_to_column('individualID')
    head(formatted_pcs)

    # make sure pheno and meth are in the same order
    head(pheno)
    ordered_pheno <- pheno %>%
        rownames_to_column('rowname') %>%
        left_join(formatted_pcs) %>%
        column_to_rownames('individualID') 
    head(ordered_pheno)
    dim(ordered_pheno)


    # create a trait column
    trait
    ordered_pheno['trait'] <- ordered_pheno[trait] 
    head(ordered_pheno)

    # order the date to connect them
    ordered_pheno <- ordered_pheno %>%
        filter(!is.na(trait)) 
    dim(ordered_pheno)

    ordered_samples <- intersect(rownames(ordered_pheno), colnames(meth))
    head(ordered_samples)

    # move gene id to rownames
    # and normalize
    head(meth)
    ordered_meth <- meth[,ordered_samples] %>%
        na.omit() %>%
        apply(1,RankNorm) %>%
        t() %>%
        as.data.frame()
    head(ordered_meth)

    ordered_pheno <- ordered_pheno[ordered_samples,]

    dim(ordered_meth)
    dim(ordered_pheno)

    # check to make sure they're ordered identically
    if(!identical(colnames(ordered_meth), rownames(ordered_pheno))){
      stop("ERROR: colnames mvals don't match targets")
    }

        
    # ---- DMP Analysis with limma ----
    # create model

    head(ordered_pheno)

    # replace missing pmi with mean
    ordered_pheno <- ordered_pheno %>%
        mutate(pmi = ifelse(is.na(pmi), mean(pmi, na.rm=TRUE), pmi)) %>%
        mutate(age_death = ifelse(is.na(age_death), mean(age_death, na.rm=TRUE), age_death)) 

        


    # set up the model
    model_content <-  paste0("~ trait + age_death + sex + pmi + batch + study")

    # include global pcs or not?
    if (include_global){
        model_content <-  paste0(paste(c(model_content, pc_cols), collapse=' + '))
    } 

    # add cell type proportions
    model_content <- paste0(paste(c(model_content, c("neuron", "endo", "astro")), collapse=" + "))
    model_content

    if (test_bulk){
        model_content <- paste0(c("neuron", "oligo_opc", "endo", "astro", "age_death + sex + pmi + batch + study", pc_cols), collapse=" + ")
        model_content <- paste0("~ ", model_content)
    }
    model_content

    print(model_content)


    # design model 
    design <- model.matrix(as.formula(model_content), data = ordered_pheno)

    head(design)
    dim(ordered_meth)
    dim(design)
        
    # run linear model
    lmfit <- lmFit(ordered_meth, design)

    # remove number of unprotected svs from degrees of freedom in lmfit
    # do not need to
    #remove <- 0
    #lmfit$df.residual <- lmfit$df.residual - remove
    
    # use empirical bayes to stabilize estimates
    lm <- eBayes(lmfit)
        
    # get number of probes that are significant with FDR <0.05
    if (!test_bulk){
        head(lm)
        results <- topTable(
          lm,
          coef = 'trait', 
          number = Inf,
          adjust.method = 'none'
        )
        head(results)
    }

    res <- lapply(c('neuron', 'oligo_opc', 'endo', 'astro'), get_res)
    results <- do.call(rbind, res)
    head(results)

    # add p-value adjustment here
    results <- results %>%
        mutate(bh_pval = p.adjust(P.Value, method='BH'),
               bf_pval = p.adjust(P.Value, method='bonferroni')
        )
        
    # filter results
    results_sub <- results[results$bh_pval < 0.05, ]

    dim(results)
    dim(results_sub)

    print(paste("Number of sig results:", nrow(results_sub), "out of", nrow(results)))

    # save the filtered results
    filename <- paste(trait, cell_type, region_type, summary_type, "filtered_limma_results.rds", sep='_')
    savefile <- paste0(savedir, cleaned_str, filename)
    savefile
    saveRDS(results_sub, savefile)

    1
}

load_meth <- function(region_type, summary_type, cleaned_str, cell_type){
    # load the methylation data

    # set up the name of the data file
    datafile <- NA
    if (summary_type == 'cpgs'){
        datafile <- paste0(datadir, precleaned_str, cleaned_str, cell_type, "_lifted_hg38_mvals.rds")
    } else{
        datafile <- paste0(datadir, cleaned_str, region_type, "_", cell_type, "_", summary_type, ".rds")
    }
    datafile

    # read in the data file
    meth <- readRDS(datafile)
    head(meth)[1:10]
    dim(meth)


    # format the data so they're uniform across summary types and cell types
    if (summary_type == 'avgs' & region_type != 'astro_enh'){
        meth <- meth %>%
            select(-duration_secs)
    }
    if (summary_type == 'pcs'){
        # remove columns that are not pcs
        head(meth)
        remove_cols <- c('num_cpgs', 'gv_dim', 'mp_dim', 'duration_secs', 'error', 'percent_variance', 'percent_variance_explained')

        keep_cols <- setdiff(colnames(meth), remove_cols)
        length(keep_cols)
        dim(meth)

        sub_meth <- meth[,keep_cols] %>%
            na.omit()
        head(sub_meth)

        meth <- sub_meth

    }

    meth
}

check_results <- function(){
    # check the results and get some metrics

    # define parameter choices for all runs
    region_types <- c('gene_body', 'preTSS', 'full_gene', 'promoters')
    summary_types <- c('pcs', 'avgs', 'cpgs')
    traits <- c('ceradsc', 'cat_braaksc', 'new_apoe', 'new_diag')
    cleaned_strs <- c("cleaned_global_")
    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')

    runs <- expand.grid(region_type=region_types,
                        #summary_type=summary_types,
                        trait=traits,
                        cleaned_str=cleaned_strs,
                        cell_type=cell_types
    )
    head(runs)

    # load gene map with symbols
    datafile <- paste0("/path/to/table/gene_biomart_table.csv")
    gene_map <- read_csv(datafile) %>%
        rename(gene_no_vers=gene_id,
               gene_id=ensembl_gene_id_version
        )
    head(gene_map)


    row = runs[1,] # for testing interactively
    st = summary_types[1] # for testing interactively
    format_results <- function(row, st){
        # get parameters for this run
        region_type <- row[['region_type']]
        trait <- row[['trait']]
        cleaned_str <- row[['cleaned_str']]
        cell_type <- row[['cell_type']]

        print(paste("Processing", st, cell_type, region_type, trait, cleaned_str, Sys.time()))

        # load the dmp results
        filename <- paste(trait, cell_type, region_type, st, "unfiltered_limma_results.rds", sep='_')
        datafile <- paste0(savedir, cleaned_str, filename)
        res <- readRDS(datafile)
        head(res)

        # separate PCs from genes
        # separate enhancer types
       formatted_res <- res %>%
            rownames_to_column('feature') %>%
            separate(feature, c('gene_id', 'pc'), sep='-PC', remove=FALSE) %>%
            separate(gene_id, c('gene_id', 'enhancer_type'), sep='-') 
        head(formatted_res)

        # -- add gene symbols
        res_symbols <- formatted_res %>%
            left_join(gene_map)
        head(res_symbols)

        filename <- paste(trait, cell_type, region_type, st, "unfiltered_formatted_limma_results.rds", sep='_')
        datafile <- paste0(savedir, cleaned_str, filename)
        datafile
        saveRDS(res_symbols, datafile)

        1
    }

    apply(runs, 1, format_results, summary_types[1])
    apply(runs, 1, format_results, summary_types[2])


    # do cpgs overlap the regions that were detected using PCs?
    runs
    row <- runs[1,]
    map_cpgs <- function(row){
        summary_type <- summary_types[3]
        region_type <- row[['region_type']]
        trait <- row[['trait']]
        cleaned_str <- row[['cleaned_str']]
        cell_type <- row[['cell_type']]

        print(paste("Processing", cell_type, region_type, trait, cleaned_str, Sys.time()))

        filename <- paste(trait, cell_type, region_type, summary_type, "unfiltered_limma_results.rds", sep='_')
        datafile <- paste0(savedir, cleaned_str, filename)
        datafile
        cpg_res <- readRDS(datafile)
        head(cpg_res)

        # need to get the cpg position annotations
        mapfile <- paste0(datadir, region_type, "_", cell_type, "_gene_cpg_map.rds")
        cpg_map <- readRDS(mapfile)
        head(cpg_map)

        # need a cpgs position map
        datafile <- paste0(datadir, cell_type, "_cpg_position_map.rds")
        datafile
        cpg_position_map <- readRDS(datafile)
        head(cpg_position_map)

        full_map <- cpg_position_map %>%
            select(probe=cpg_gr37, cpg_id=cpg) %>%
            left_join(cpg_map) %>%
            filter(!is.na(gene_id))
        head(full_map)

        
        # find out which genes these cpgs map to
        head(gene_map)
        head(cpg_res)
        dim(cpg_res)
        formatted_cpg <- cpg_res %>%
            rownames_to_column('cpg_id') %>%
            left_join(full_map) %>%
            #filter(!is.na(gene_id)) %>%
            mutate(trait=trait,
                   region_type=region_type,
                   summary_type=summary_type,
                   cleaned_str=cleaned_str,
                   cell_type=cell_type
            ) %>%
            separate(gene_id, c('gene_id', 'enhancer_type'), sep='-', fill='right') 
            #separate(gene_id, c('gene_id', 'vers'), sep='\\.') %>%
            #left_join(gene_map)
        head(formatted_cpg)


        savefile <- paste(trait, cell_type, region_type, summary_type, "mapped_unfiltered_limma_results.rds", sep='_')
        savepath <- paste0(savedir, cleaned_str, savefile)
        savepath
        saveRDS(formatted_cpg, savepath)

        1
    }

    cpg_res <- apply(runs, 1, map_cpgs)


    ## -- map full array results to region types
    runs
    region_type <- region_types[4]
    row <- runs[1,]
    map_full_array_cpgs <- function(row){
        summary_type <- summary_types[3]
        region_type <- row[['region_type']]
        trait <- row[['trait']]
        cleaned_str <- row[['cleaned_str']]
        cell_type <- row[['cell_type']]

        print(paste("Processing", cell_type, region_type, trait, cleaned_str, Sys.time()))

        # read in the full array results
        filename <- paste(trait, cell_type, "full_array", summary_type, "unfiltered_limma_results.rds", sep='_')
        datafile <- paste0(savedir, cleaned_str, filename)
        datafile
        cpg_res <- readRDS(datafile)
        head(cpg_res)

        # need to get the cpg position annotations
        mapdir <- paste0("../../output/03_summarised_genes/01_summarised_genes/restimate_proportions/", cov_str, "/")
        mapfile <- paste0(mapdir, region_type, "_", cell_type, "_gene_cpg_map.rds")
        cpg_map <- readRDS(mapfile)
        head(cpg_map)

        # need a cpgs position map
        datafile <- paste0(mapdir, cell_type, "_cpg_position_map.rds")
        datafile
        cpg_position_map <- readRDS(datafile)
        head(cpg_position_map)

        full_map <- cpg_position_map %>%
            select(probe=group_name, cpg_id=cpg) %>%
            left_join(cpg_map) %>%
            filter(!is.na(gene_id))
        head(full_map)

        
        # find out which genes these cpgs map to
        head(cpg_res)
        dim(cpg_res)
        formatted_cpg <- cpg_res %>%
            rownames_to_column('cpg_id') %>%
            left_join(full_map) %>%
            #filter(!is.na(gene_id)) %>%
            mutate(trait=trait,
                   region_type=region_type,
                   summary_type=summary_type,
                   cleaned_str=cleaned_str,
                   cell_type=cell_type
            ) %>%
            separate(gene_id, c('gene_id', 'enhancer_type'), sep='-') %>%
            #separate(gene_id, c('gene_id', 'vers'), sep='\\.') %>%
            left_join(gene_map)
        head(formatted_cpg)


        savefile <- paste(trait, cell_type, region_type, summary_type, "full_array_mapped_unfiltered_limma_results.rds", sep='_')
        savepath <- paste0(savedir, cleaned_str, savefile)
        savepath
        saveRDS(formatted_cpg, savepath)

        1
    }

    # do not need to run
    #cpg_res <- apply(runs, 1, map_full_array_cpgs)


    1
}

match_cpgs <- function(meth, region_type, cleaned_str, cell_type){
    # match cpgs with those used for pcs

    # need a cpgs position map
    datafile <- paste0(mapdir, cell_type, "_cpg_position_map.rds")
    datafile
    cpg_position_map <- readRDS(datafile)
    head(cpg_position_map)

    if (region_type == 'full_array'){
        # just switch out cpg names
        head(meth)

        # subset to cpgs in meth data
        sub_map <- cpg_position_map %>%
            filter(group_name %in% rownames(meth))
        head(sub_map)
        dim(sub_map)
        dim(meth)

        # match ordering
        ordered_meth <- meth[sub_map$group_name,]
        head(ordered_meth)

        if (!identical(rownames(ordered_meth), sub_map$group_name)) stop("order doesn't match")

        # assign new names
        head(sub_map)
        rownames(ordered_meth) <- sub_map$cpg
        head(ordered_meth)

        return(ordered_meth)
    }

    # need to get the cpg position annotations
    mapdir <- paste0("../../output/episcore/03_summarised_genes/", cov_str)
    mapfile <- paste0(mapdir, region_type, "_", cell_type, "_gene_cpg_map.rds")
    mapfile
    cpg_map <- readRDS(mapfile)
    head(cpg_map)
    

    head(meth)
    dim(meth)

    head(cpg_map)
    head(cpg_position_map)

    # keep only cpgs that map to genes
    # and only genes that are used for pcs
    length(intersect(cpg_map$cpg_id, cpg_position_map$cpg))
    dim(cpg_map)
    #head(pc_genes)

    gene_probe_map <- cpg_position_map %>%
        rename(cpg_id=cpg) %>%
        left_join(cpg_map) %>%
        filter(!is.na(gene_id)) 
        #filter(gene_id %in% pc_genes$gene_id)
    head(gene_probe_map)
    dim(gene_probe_map)

    dim(meth)

    gene_meth <- meth[rownames(meth) %in% gene_probe_map$group_name,]
    dim(gene_meth)
    head(gene_meth)

    gene_meth
}

compute_global_pcs <- function(meth, region_type, cleaned_str, trait, cell_type){
    # compute global pcs

    savefile <- paste0(savedir, cleaned_str, trait, "_", region_type, "_", cell_type, "_global_pcs.rds")
    if (file.exists(savefile)){
        print("global pcs already saved")
        return(1)
    }


    head(meth)

    # make sure no zero variance cpgs
    var_mvals <- meth[apply(meth, 1, var) != 0,] %>%
        t()
    dim(meth)
    dim(var_mvals)

    # center and scale
    scaled_mvals <- scale(var_mvals)
    head(scaled_mvals)[1:10]
    dim(scaled_mvals)

    # estimate dimension
    # rows = features, cols=samples
    ndim <- EstDimRMT(t(scaled_mvals))$dim
    ndim
    
    # run pca
    pc_res <- prcomp(t(scaled_mvals))
    names(pc_res)
    
    pcs = pc_res$rotation[,1:ndim,drop=F]
    head(pcs)

    savefile <- paste0(savedir, cleaned_str, trait, "_", region_type, "_", cell_type, "_global_pcs.rds")
    savefile
    saveRDS(pcs, savefile)

    1
}

## --- plotting functions -- ##
plot_dmp_counts_barplot <- function(combined_df, correction_type, summary_type_colors, plotting_dir, full_array){
    head(combined_df)

    # select the method that was run
    # select the path that matches
    full_array_str <- "full_array"
    if (!full_array) full_array_str <- "matched_array"
    sub_dir <- paste0(plotting_dir, correction_type, "_", full_array_str, "/")
    dir.create(sub_dir, showWarnings=FALSE)

    # set the p-value correction type
    combined_df$adj_pval <- combined_df$bh_pval
    if (correction_type == 'bf') combined_df$adj_pval <- combined_df$bf_pval

    # change feature level of cpgs to just the probes
    # rather than gene-probe pairs
    head(combined_df)
    formatted_cpgs <- combined_df %>%
        filter(summary_type == 'CpGs') %>%
        separate(feature_id, c('gene_id', 'feature'), sep='-') %>%
        select(-gene_id) %>%
        arrange(feature) %>%
        unique()
    head(formatted_cpgs)

    # combined these formatted cpgs to the other summary types
    head(combined_df)
    formatted_df <- combined_df %>%
        filter(summary_type != 'CpGs') %>%
        full_join(formatted_cpgs)
    head(formatted_df)

    # filter for using the full array or not
    if (full_array){
        formatted_df <- formatted_df %>%
            filter(full_array)
    } else{
        formatted_df <- formatted_df %>%
            filter(!full_array)
    }

    # filter fo significant p-vals
    # count the number of sig reults for each group
    head(formatted_df)
    sig_counts <- formatted_df %>%
        filter(adj_pval < 0.05) %>%
        select(feature_id, trait, summary_type, region_type, cell_type) %>%
        unique() %>%
        count(region_type, trait, summary_type, cell_type) %>%
        rename(sig_dmp=n)
    head(sig_counts)


    # count the total number of input features for each group
    head(formatted_df)
    total_counts <- formatted_df %>%
        select(feature_id, trait, summary_type, region_type, cell_type) %>%
        unique() %>%
        count(region_type, trait, summary_type, cell_type) 
    head(total_counts)

    # create a wide data frame of counts
    wide_counts <- total_counts %>%
        spread(cell_type, n) %>%
        arrange(summary_type)
    wide_counts


    # compute the percentage of sig results for each group
    head(sig_counts)
    sig_counts <- sig_counts %>%
        left_join(total_counts) %>%
        mutate(percent = sig_dmp / n * 100) %>%
        mutate(label = paste0(sig_dmp))
    head(sig_counts)

    # change the region type for cpgs
    # set the position for text label
    # set color for summary types
    formatted_counts <- sig_counts %>%
        unique() %>%
        ungroup() %>%
        mutate(max_dmp = max(sig_dmp)) %>%
        mutate(ypos = ifelse(sig_dmp < max_dmp-(max_dmp/2),
                             sig_dmp,
                             sig_dmp-(max_dmp*(1/2)))) %>%
        mutate(color = ifelse(sig_dmp < max_dmp-(max_dmp/2),
                              summary_type, 'white'))
    head(formatted_counts)

    compare_feature_proportions <- function(){
        # compare proportion of features that were found to be significant
        # out of the total number of features
        head(formatted_df)

        cmp_df <- formatted_df %>%
            mutate(dm = adj_pval < 0.05) %>%
            select(feature_id, summary_type, region_type, cell_type, adj_pval, trait, dm)
        head(cmp_df)

        ad_counts <- cmp_df

        runs <- expand.grid(cell_type = unique(ad_counts$cell_type),
                       region_type = unique(ad_counts$region_type),
                       trait = unique(ad_counts$trait),
                       st1 = unique(ad_counts$summary_type),
                       st2 = unique(ad_counts$summary_type)
        ) %>%
            filter(
                   st1 == 'pcs'
            ) %>%
            mutate(rn = row_number()) %>%
            filter(region_type == 'full_gene')
        head(runs)
        dim(runs)

        row <- runs[1,]
        test_feature_props <- function(row){
            # define the parameters for this run
            ct <- row[['cell_type']] 
            rt <- row[['region_type']] %>% as.character()
            tr <- row[['trait']] %>% as.character()
            st1 <- row[['st1']] %>% as.character()
            st2 <- row[['st2']] %>% as.character()
            rn <- row[['rn']] %>% as.character()

            # filter for the parameters
            row_df <- ad_counts %>%
                filter(region_type == rt,
                       trait == tr,
                       cell_type == ct, 
                       summary_type %in% c(st1, st2)
                ) %>%
                mutate(summary_type = as.character(summary_type)) 
            head(row_df)

            # checks to avoid errors
            if (length(unique(as.character(row_df$summary_type))) < 2) return(NA)
            if (length(unique(as.character(row_df$dm))) < 2) return(NA)

            print(paste(st1, st2, ct, rt, tr, rn))

            # reate contingency table
            tab <- table(row_df$dm, row_df$summary_type)
            tab

            # run a fishers test
            res <- fisher.test(tab)
            res
            names(res)

            # return the results
            data.frame(region_type=rt,
                       trait=tr,
                       cell_type=ct,
                       st1 = st1,
                       st2 = st2,
                       pval=res$p.value,
                       ndm1 = tab['FALSE',st1],
                       ndm2 = tab['FALSE',st2],
                       dm1 = tab['TRUE',st1],
                       dm2 = tab['TRUE',st2]
            )
        }
        # run for all combinations of parameters
        props <- apply(runs , 1, test_feature_props)
        props_df <- do.call(rbind, props) %>%
            na.omit()
        head(props_df)

        # correct for multiple tests
        props_corrected <- props_df %>%
            group_by(trait, st1, st2, region_type) %>%
            mutate(adj_pval = p.adjust(pval, method='bonferroni')) %>%
            mutate(prop1 = dm1 / (dm1+ndm1),
                   prop2 = dm2 / (dm2+ndm2)
            ) %>%
            mutate(sig = ifelse(adj_pval < 0.1, '.', ''),
                   sig = ifelse(adj_pval < 0.05, '*', sig),
                   sig = ifelse(adj_pval < 0.01, '**', sig),
                   sig = ifelse(adj_pval < 0.001, '***', sig),
            ) %>%
            arrange(cell_type)
        head(props_corrected)

        table(props_corrected$pval < 0.05)
        table(props_corrected$adj_pval < 0.05)

        # format the data frame of restults
        head(props_corrected)
        formatted_props <- props_corrected %>%
            select(region_type, trait, cell_type, st1, st2, adj_pval, prop1, prop2, sig) %>%
            gather('prop_type', 'prop_dm', prop1:prop2) %>%
            mutate(summary_type = ifelse(prop_type == 'prop1', st1, st2),
                   new_sig = ifelse(summary_type == 'PCs', '', sig)) %>%
            ungroup() %>%
            select(-prop_type, -st1, -st2, -adj_pval, -sig) %>%
            unique() %>%
            left_join(cell_type_map) %>%
            left_join(summary_type_map) 
        head(formatted_props)

        formatted_sig <- formatted_props %>%
            group_by(region_type, cell_type) %>%
            mutate(max_prop=max(prop_dm)) %>%
            select(-prop_dm) %>%
            unique()
        head(formatted_sig)

        # plot the results as a bar plot
        head(formatted_props)
        p <- formatted_props %>%
            ggplot() +
            geom_bar(aes(x=formatted_ct, y=prop_dm, fill=formatted_st),
                     stat='identity', 
                     position='dodge'
            ) +
            geom_text(data=formatted_props, aes(x=formatted_ct, y=prop_dm+0.001, label=new_sig, color=formatted_st),
                      size=5, position=position_dodge(width=0.9)) +
            scale_fill_manual(values=summary_type_colors, name="Summary type") +
            scale_color_manual(values=summary_type_colors, name="Summary type") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1),
                  text = element_text(size=20),
                  legend.position='top',
                  legend.justification='right'
            ) +
            xlab("Cell type") +
            ylab("Proportion of DM features") 

        # save the plot and the rds object
        trait <- 'ceradsc'
        (savefile <- paste0(sub_dir, trait, "_feature_proportions_comparison.png"))
        ggsave(p, file=savefile, width=7, height=7)

        (savefile <- paste0(sub_dir, trait, "_feature_proportions_comparison.rds"))
        saveRDS(p, savefile)

        # plot the total counts of features also
        head(props_corrected)
        total_counts <- props_corrected %>%
            ungroup() %>%
            gather('ndm_type', 'ndm', ndm1:ndm2) %>%
            gather('dm_type', 'dm', dm1:dm2) %>%
            gather('st_type', 'st', st1:st2) %>%
            mutate(ndm_type = str_remove(ndm_type, "ndm"),
                   dm_type = str_remove(dm_type, 'dm'),
                   st_type = str_remove(st_type, 'st')
            ) %>%
            filter(ndm_type == dm_type, dm_type == st_type) %>%
            mutate(total_genes = ndm + dm) %>%
            select(cell_type, st, total_genes) %>%
            unique() %>%
            rename(summary_type = st) %>%
            group_by(summary_type) %>%
            summarize(total_genes = median(total_genes))
        head(total_counts)

        p <- ggplot(total_counts, aes(x=summary_type, y=total_genes)) +
            geom_bar(aes(fill=summary_type),
                     stat='identity', 
                     position='dodge'
            ) +
            geom_label(aes(label=total_genes, color=summary_type),
                      size=5, position=position_dodge(width=0.9)) +
            #facet_wrap(region_type ~ ., nrow=2) +
            scale_fill_manual(values=summary_type_colors, name="Summary type") +
            scale_color_manual(values=summary_type_colors, name="Summary type") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 50, hjust=1, vjust=1),
                  text = element_text(size=20), 
                  legend.position='none'
            ) +
            xlab("Cell type") +
            ylab("Total features") +
            ggtitle(trait)
        savefile <- paste0(sub_dir, trait, "_feature_total_comparison.png")
        ggsave(p, file=savefile, width=7, height=7)
    }


    # plot the counts of significant DMPs with traits
    with_cpgs=TRUE
    plot_full_gene_counts_with_traits <- function(with_cpgs){
        head(formatted_counts)

        filtered_counts <- formatted_counts %>%
            filter(summary_type != 'CpGs') %>%
            filter(region_type == 'Full gene') 
        if (with_cpgs){
            filtered_counts <- formatted_counts %>%
                filter(region_type == 'Full gene') 
        }
        head(filtered_counts)

        p <- filtered_counts %>%
            mutate(sig_genes = sig_dmp) %>%
            ggplot(aes(x=cell_type, y=sig_genes, fill=summary_type)) +
            geom_bar(stat='identity', position='dodge') +
            geom_text(aes(label=label,
                          y=ypos,
                          color=color
                          ),
                      hjust=0.5,
                      vjust=-0.5,
                      show.legend=FALSE,
                      size=5,
                      position=position_dodge(width = .9)
                      ) +
            facet_grid(trait ~ ., scales='fixed') +
            scale_color_manual(values=summary_type_colors, breaks=summary_types) +
            scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
            theme_bw() +
            theme(
                  strip.background=element_rect(fill='white', color='black'),
                  panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  text=element_text(size=20, family='Open Sans'),
                  legend.title=element_text(size=20, face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
            ) +
            #coord_cartesian(ylim=c(0,270)) +
            ggtitle(paste(correction_type, "sig genes")) +
            xlab("Cell types") +
            ylab("Significant DM probes")

            main_filename = "_traits_fullGene_dm_DMP_counts_barplot.png"
            # not including cpgs
            savefile <- paste0(sub_dir, correction_type, 
                               main_filename)
            if (full_array){
                savefile <- paste0(sub_dir, correction_type, 
                                   "_fullArray", 
                                   main_filename)
            }
            if (with_cpgs){
                # inlucde cpgs
                savefile <- paste0(sub_dir, correction_type, 
                                   "_withCpgs",
                                   main_filename)
                if (full_array){
                    savefile <- paste0(sub_dir, correction_type, 
                                       "_withCpgs",
                                       "_fullArray", 
                                       main_filename)
                }
            }
            savefile

            ggsave(savefile, plot=p, width=12, height=7, units='in')
        1
    }
    
    plot_full_gene_counts_with_traits(with_cpgs=FALSE)
    plot_full_gene_counts_with_traits(with_cpgs=TRUE)


    # plot the counts of significant DMPs
    head(formatted_counts)
    p <- formatted_counts %>%
        filter(summary_type %in% c('Avgs', 'PCs')) %>%
        ggplot(aes(x=summary_type, y=sig_dmp, fill=summary_type)) +
        geom_bar(stat='identity', position='dodge') +
        geom_text(aes(label=label,
                      y=ypos,
                      color=color
                      ),
                  hjust=0.5,
                  vjust=-0.5,
                  show.legend=FALSE,
                  size=10,
                  ) +
        facet_grid(region_type ~ cell_type, scales='fixed') +
        scale_color_manual(values=summary_type_colors, breaks=summary_types) +
        scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
        theme_bw() +
        theme(
              text=element_text(size=20, family='Open Sans'),
              legend.title=element_text(size=20, face='bold'),
              plot.title=element_text(hjust=0.5),
              axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
        ) +
        ggtitle(paste(trait, correction_type, "Sig DMPs")) +
        xlab("Summary types") +
        ylab("Number of significant DM regions")


    savefile <- paste0(sub_dir, trait, "_", correction_type, 
                       "_dmp_counts_barplot.png")
    if (full_array){
        savefile <- paste0(sub_dir, trait, "_", correction_type, 
                           "_full_array_dmp_counts_barplot.png")
    }
    savefile

    ggsave(savefile, plot=p, width=18, height=13, units='in')


    # show only partial region types
    head(formatted_counts)
    trait <- "ceradsc"
    p <- formatted_counts %>%
        #filter(region_type %in% c('Full genes', 'Promoters', 'CpGi')) %>%
        ggplot(aes(x=summary_type, y=sig_dmp, fill=summary_type)) +
        geom_bar(stat='identity', position='dodge') +
        geom_text(aes(label=label,
                      y=ypos,
                      color=color
                      ),
                  hjust=0.5,
                  vjust=-0.5,
                  show.legend=FALSE,
                  size=10,
                  ) +
        facet_grid(region_type ~ cell_type, scales='fixed') +
        scale_color_manual(values=summary_type_colors, breaks=summary_types) +
        scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
        theme_bw() +
        theme(
              text=element_text(size=20, family='Open Sans'),
              legend.title=element_text(size=20, face='bold'),
              plot.title=element_text(hjust=0.5),
              axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
        ) +
        ggtitle(paste(trait, correction_type, "Sig DMPs")) +
        xlab("Summary types") +
        ylab("Number of significant DM regions")


    savefile <- paste0(sub_dir, trait, "_", correction_type, 
                       "_sub_dmp_counts_barplot.png")
    if (full_array){
        savefile <- paste0(sub_dir, trait, "_", correction_type, 
                           "_full_array_sub_dmp_counts_barplot.png")
    }
    savefile

    ggsave(savefile, plot=p, width=18, height=10, units='in')



    # plot the percentages of the counts
    head(formatted_counts)
    tr="ceradsc"
    ymax = max(formatted_counts$percent)
    p <- formatted_counts %>%
        filter(region_type == 'Full gene',
               trait == tr) %>%
        filter(full_array = TRUE) %>%
        mutate(label=paste0(round(percent, 2))) %>%
        ggplot(aes(x=summary_type, y=percent, fill=summary_type)) +
        geom_bar(stat='identity', position='dodge') +
        geom_text(aes(label=label,
                      y=percent,
                      color=summary_type
                      ),
                  color='black',
                  hjust=0.5,
                  vjust=-0.5,
                  show.legend=FALSE,
                  size=7,
                  position=position_dodge(width=1)
                  ) +
        facet_grid(region_type ~ cell_type, scales='fixed') +
        #scale_color_manual(values=summary_type_colors, breaks=summary_types) +
        scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
        scale_color_manual(values=summary_type_colors, breaks=summary_types) +
        theme_bw() +
        theme(
              text=element_text(size=20, family='Open Sans'),
              legend.title=element_text(size=20, face='bold'),
              plot.title=element_text(hjust=0.5),
              axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
        ) +
        ggtitle(paste(trait, correction_type, "Sig DMPs")) +
        xlab("Summary types") +
        ylab("Percentage of significant DM regions") +
        coord_cartesian(ylim=c(0,ymax+0.015))

    savefile <- paste0(sub_dir, trait, "_", correction_type, 
                       "_dmp_percentage_barplot.png")
    if (full_array){
        savefile <- paste0(sub_dir, trait, "_", correction_type, 
                           "_full_array_dmp_percentage_barplot.png")
    }
    savefile

    ggsave(savefile, plot=p, width=11, height=8, units='in')


    # plot percentage with region types
    head(formatted_counts)
    head(props_corrected)
    tr="ceradsc"
    ct = "Bulk"
    filtered_counts <- formatted_counts %>%
        mutate(lesser_st = as.character(summary_type)) %>%
        left_join(props_corrected) %>%
        filter(greater_st == "PCs" | summary_type == "PCs"
        )  %>%
        filter(cell_type == ct,
               trait == tr) %>%
        mutate(sig_label = ifelse(pval < 0.1, ".", NA),
               sig_label = ifelse(pval < 0.05, "*", sig_label),
               sig_label = ifelse(pval < 0.001, "**", sig_label),
        )
    head(filtered_counts)


    ymax = max(filtered_counts$percent)
    p <- filtered_counts %>%
        filter(full_array = TRUE) %>%
        mutate(label=paste0(round(percent, 2))) %>%
        ggplot(aes(x=region_type, y=percent, fill=summary_type)) +
        geom_bar(stat='identity', position='dodge') +
        geom_text(aes(label=sig_label,
                      y=percent,
                      color=summary_type
                      ),
                  #color='black',
                  hjust=0.5,
                  vjust=-0.3,
                  #show.legend=FALSE,
                  size=6,
                  position=position_dodge(width=1)
                  ) +
        facet_grid(trait ~ ., scales='fixed') +
        scale_fill_manual(values=summary_type_colors, breaks=names(summary_type_colors)) +
        scale_color_manual(values=summary_type_colors, breaks=names(summary_type_colors)) +
        theme_bw() +
        theme(
              text=element_text(size=20, family='Open Sans'),
              legend.title=element_text(size=20, face='bold'),
              plot.title=element_text(hjust=0.5),
              axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
              axis.text.x=element_text(angle=30, vjust=1, hjust=1)
        ) +
        labs(fill=paste("Summary type"), color=paste("Summary type")) +
        ggtitle(paste(trait, correction_type, "Sig DMPs")) +
        xlab("Summary types") +
        ylab("Percentage of significant DM regions") +
        coord_cartesian(ylim=c(0,ymax+0.015))

    savefile <- paste0(sub_dir, trait, "_", correction_type, 
                       "_allRegions",
                       "_dmp_percentage_barplot.png")
    if (full_array){
        savefile <- paste0(sub_dir, trait, "_", correction_type, 
                           "_full_array_dmp_percentage_barplot.png")
    }
    savefile

    ggsave(savefile, plot=p, width=11, height=7, units='in')



    ## -- plot the number of total features per summary type
    head(sig_counts)
    counts <- sig_counts %>%
        group_by(region_type) %>%
        rename(total_probes=n) %>%
        mutate(max_dmp=max(total_probes)) %>%
        mutate(ypos=ifelse(total_probes < max_dmp-(max_dmp/2),
                           total_probes,
                           total_probes-max_dmp*(1/2))) %>%
        mutate(color=ifelse(total_probes < max_dmp-(max_dmp/2),
                            as.character(summary_type),
                            'black'))
    head(counts)

    summary_type_colors['black'] <- 'black'

    p <- counts %>%
        ggplot(aes(x=summary_type, y=total_probes, fill=summary_type)) +
        geom_bar(stat='identity', position='dodge') +
        geom_text(aes(label=total_probes,
                      y=ypos,
                      color=color
                      ),
                  hjust=0.5,
                  vjust=-0.5,
                  show.legend=FALSE,
                  size=10
                  #position=position_dodge(width=1)
                  ) +
        facet_grid(region_type ~ cell_type, scales='free_y') +
        scale_color_manual(values=summary_type_colors, breaks=summary_types) +
        scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
        theme_bw() +
        theme(
              text=element_text(size=20, family='Open Sans'),
              legend.title=element_text(size=20, face='bold'),
              plot.title=element_text(hjust=0.5),
              axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
        ) +
        ggtitle(paste(trait, correction_type, "Sig DMPs")) +
        xlab("Summary types") +
        ylab("Number of total features")

    savefile <- paste0(sub_dir, trait, "_", correction_type, 
                       "_feature_count_barplot.png")
    if (full_array){
        savefile <- paste0(sub_dir, trait, "_", correction_type, 
                           "_full_array_feature_count_barplot.png")
    }
    savefile

    ggsave(savefile, plot=p, width=18, height=13, units='in')


    1
}


plot_gene_counts_barplot <- function(combined_df, correction_type, summary_type_colors, plotting_dir, full_array){
    # plot results on a gene-level
    head(combined_df)

    print(paste("correction type", correction_type))

    full_array_str <- "full_array"
    if (!full_array) full_array_str <- "matched_array"

    sub_dir <- paste0(plotting_dir, correction_type, "_", full_array_str, "/")
    dir.create(sub_dir, showWarnings=FALSE)

    # set the p-value correction type
    combined_df$adj_pval <- combined_df$bh_pval
    if (correction_type == 'bf') combined_df$adj_pval <- combined_df$bf_pval

    # filter for using the full array or not
    if (full_array){
        combined_df <- combined_df %>%
            filter(full_array)
    } else{
        combined_df <- combined_df %>%
            filter(!full_array)
    }

    # change feature level of cpgs to just the probes
    # rather than gene-probe pairs
    formatted_df <- combined_df %>%
        mutate(full_array = TRUE)
    head(formatted_df)

    check_sig <- function(){
        # get number of cpgs per pc
        list.files(datadir)
        filename <- paste(region_type, cell_type, "gene_cpg_map.rds", sep='_')
        datafile <- paste0(datadir, filename)
        cpg_map <- readRDS(datafile)
        head(cpg_map)

        # get the number of cpgs per gene
        cpg_counts <- cpg_map %>%
            count(gene_id, name="num_cpgs")
        head(cpg_counts)

        # filter for significant results
        sig_df <- formatted_df %>%
            filter(adj_pval < 0.05) %>%
            arrange(adj_pval) %>%
            left_join(cpg_counts) %>%
            select(gene_id, num_cpgs) %>%
            unique()
        table(sig_df$feature)
        head(sig_df)

        # filter for non-isgnificant results
        nsig_df <- formatted_df %>%
            filter(adj_pval > 0.05) %>%
            left_join(cpg_counts) %>%
            select(gene_id, num_cpgs) %>%
            unique() %>%
            filter(!gene_id %in% sig_df$gene_id)
        head(nsig_df)

        summary(sig_df$num_cpgs)
        summary(nsig_df$num_cpgs)

    }



    # check genes that were significant in INTACT results
    check_intact_genes <- function(){
        head(formatted_df)
        intact_genes <- c('MS4A6E', 'MS4A6A', 'SPACDR', 'APOC1', 'RNF5', 'APOE', 'APOC4',
                          'TSBP1', 'ERCC2', 'RELB', 'PSMG3', 'DEF6', 'MS4A4A', 'RNF111',
                          'ILRUN', 'SLC39A13', 'SEPTIN4', 'ADGRF4', 'MAMSTR', 'STAG3'
        )

        # load gene symbols or create the file
        savefile <- paste0(savedir, "gene_symbols.rds")
        if (!file.exists(savefile)){
            ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
            dm_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                     #filters = "ensembl_gene_id", 
                     #values = unique(genes_novers), 
                     mart = ensembl) 
            saveRDS(dm_symbols, savefile)
        } else{
            print("Loading existing biomart genes file")
            dm_symbols <- readRDS(savefile)
        }
        head(dm_symbols)

        # attach symbols to DM results
        head(formatted_df)
        symbols_df <- formatted_df %>%
            mutate(ensembl_gene_id = str_remove(gene_id, "\\.[0-9]+")) %>%
            left_join(dm_symbols)
        head(symbols_df)

        # filter to intact genes
        intact_dm <- symbols_df %>%
            filter(hgnc_symbol %in% intact_genes) %>%
            mutate(sig_dm = adj_pval < 0.05) %>%
            select(symbol=hgnc_symbol, sig_dm, trait, summary_type, region_type, cell_type) %>%
            group_by(symbol, trait, summary_type, region_type, cell_type) %>%
            summarize(sig_dm = any(sig_dm)) 
        head(intact_dm)
        dim(intact_dm)

        filtered_dm <- intact_dm %>%
            group_by(symbol, trait, summary_type, region_type) %>%
            mutate(sig_count = sum(as.numeric(sig_dm))) %>%
            filter(sig_count > 0) %>%
            spread(cell_type, sig_dm)
        head(filtered_dm)

    }

    # plot combined summary type upset plot
    tr = 'ceradsc'
    plot_celltype_upset_withBothSummaries <- function(tr){
        # create upset plot
        head(formatted_df)

        # formta the cell type upset input
        head(formatted_df)
        pcs_df <- formatted_df %>%
            mutate(summary_type=formatted_st,
                   cell_type=formatted_ct,
                   region_type=formatted_rt
            ) %>%
            filter(adj_pval < 0.05) %>%
            filter(summary_type == 'rPCs',
                   trait == tr) %>%
            filter(region_type == 'Full gene') %>%
            select(gene_id, cell_type) %>%
            unique() %>%
            mutate(exists = 1) %>%
            spread(cell_type, exists, fill=0) 
        head(pcs_df)

        avgs_df <- formatted_df %>%
            mutate(summary_type=formatted_st,
                   cell_type=formatted_ct,
                   region_type=formatted_rt
            ) %>%
            filter(adj_pval < 0.05) %>%
            filter(summary_type == 'Avgs',
                   trait == tr) %>%
            filter(region_type == 'Full gene') %>%
            select(gene_id, cell_type) %>%
            unique() %>%
            mutate(exists = 1) %>%
            spread(cell_type, exists, fill=0) 
        head(avgs_df)

        # get celltype names
        (cts <- colnames(avgs_df)[2:ncol(avgs_df)])

        # convert group indivator columns to boolean
        avgs_df[cts] <- avgs_df[cts] == 1
        pcs_df[cts] <- pcs_df[cts] == 1
        head(avgs_df)

        # get intersection groups for all combinations of cell types
        cts
        group_df <- avgs_df
        ints <- data.frame(expected_count=upset_data(group_df, intersect=cts)$sizes$exclusive_intersection) %>%
            arrange(-expected_count) 
            rownames_to_column('int') %>%
            mutate(cts = int %in% cts) %>%
            group_by(cts) %>%
            summarise(count=sum(expected_count))
        head(ints)

        # get intersection sizes
        group_df <- pcs_df
        get_top <- function(group_df){
            ints <- data.frame(expected_count=upset_data(group_df, intersect=cts)$sizes$exclusive_intersection) %>%
            arrange(-expected_count) %>%
            head(12)
        }
        pcs_top <- get_top(pcs_df)
        avgs_top <- get_top(avgs_df)
        pcs_top

        # take the top intersections from both summary types
        top_ints <- unique(rownames(pcs_top), rownames(avgs_top))
        top_ints

        # get the intersection group sizes
        group_df <- pcs_df
        get_counts <- function(group_df){
            ints <- data.frame(expected_count=upset_data(group_df, intersect=cts)$sizes$exclusive_intersection)
            head(ints)

            ints$intersection <- rownames(ints)
            ints <- ints[top_ints,]
            ints
        }

        avg_ints <- get_counts(avgs_df)
        pc_ints <- get_counts(pcs_df)

        head(avg_ints)
        head(pc_ints)

        # join the data together
        joined_ints <- avg_ints %>%
            as.data.frame() %>%
            rename(Avgs=expected_count) %>%
            full_join(pc_ints) %>%
            rename(rPCs=expected_count) %>%
            gather('summary_type', 'count', Avgs, rPCs) %>%
            #mutate(intersection = str_split(intersection, pattern='-')) %>%
            arrange(-count) %>%
            mutate(order = row_number()) %>%
            group_by(intersection) %>%
            mutate(order = min(order))
        head(joined_ints)

        # order the intersection groups 
        ordered_ints <- unique(joined_ints$intersection)
        joined_ints$intersection <- factor(joined_ints$intersection, levels=ordered_ints)

        # plot the results as upset plot
        (savefile <- paste0(sub_dir, tr, "_intersect_barplot.png"))
        head(joined_ints)
        p <- joined_ints %>%
            ggplot(aes(x=intersection, y=count, group=summary_type)) +
            geom_bar(aes(y=count, fill=summary_type),
                     stat='identity',
                     position='dodge'
            ) +
            scale_color_manual(values=summary_type_colors, guide='none') +
            scale_fill_manual(values=summary_type_colors, name="Summary type") +
            axis_combmatrix(sep = '-') +
            theme_bw() +
            theme(text=element_text(size=20),
                  strip.background=element_rect(fill='white', color='black'),
                  panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=1),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  legend.position=c(0.75,0.85),
            ) +
            theme_combmatrix(combmatrix.label.text=element_text(size=16)) +
            xlab("Intersection") +
            ylab("Intersection size") 
        ggsave(p, file=savefile, width=7, height=7)

        (savefile <- paste0(sub_dir, tr, "_celltype_intersect_barplot.rds"))
        saveRDS(p, savefile)

    }

    head(formatted_df)
    st = 'avgs'
    tr = 'ceradsc'
    plot_upset <- function(st, tr){
        # create upset plot
        head(formatted_df)

        # formta the cell type upset input
        upset_df <- formatted_df %>%
            filter(adj_pval < 0.05) %>%
            filter(summary_type == st,
                   trait == tr) %>%
            filter(region_type == 'full_gene') %>%
            select(gene_id, cell_type) %>%
            mutate(exists = 1) %>%
            spread(cell_type, exists, fill=0) 
        head(upset_df)

        (cts <- colnames(upset_df)[2:ncol(upset_df)])

        # convert group indivator columns to boolean
        upset_df[cts] <- upset_df[cts] == 1
        head(upset_df)


        (savefile <- paste0(sub_dir, tr, "_summarytypes_full_gene_celltype_upset.pdf"))
        pdf(savefile, width=9, height=9)
        p <- upset(upset_df, cts)
        print(p)
        dev.off()


        gene_lists <- formatted_df %>%
            filter(adj_pval < 0.05) %>%
            filter(summary_type == st,
                   trait == tr) %>%
            filter(region_type == 'full_gene') %>%
            distinct(gene_id, summary_type, cell_type) %>%
            unique() %>%
            group_by(summary_type, cell_type) %>%
            summarise(genes=list(gene_id)) 
        head(gene_lists)

        # create upset plot for the same data
        selectlists <- gene_lists$genes
        names(selectlists) <- gene_lists$cell_type
        gene_lists$cell_type
        gene_lists

        # set the colors for plotting
        cell_type_colors
        names(cell_type_colors)[2] <- "Endothelial" # update for these plots
        ordered_cts <- gene_lists %>%
            rowwise() %>%
            mutate(numgenes = length(genes)) %>%
            arrange(-numgenes) %>%
            mutate(cell_type=as.character(cell_type))
        head(ordered_cts)
        ordered_colors <- cell_type_colors[ordered_cts$cell_type]
        ordered_colors

        savefile <- paste0(sub_dir, tr, '_', st, "_full_gene_celltype_upset.pdf")
        savefile
        pdf(savefile, width=7, height=6)
        p <- upset(fromList(selectlists), 
              order.by='freq',
              nintersects=10,
              sets.bar.color=ordered_colors,
              text.scale=c(2, 2, 2, 1, 2, 2)
        )
        print(p)
        dev.off()
    }
    plot_upset(st='PCs', 'ceradsc')
    plot_upset(st='Avgs', 'ceradsc')

    head(formatted_df)
    st = 'Avgs'
    tr = 'ceradsc'
    ct = 'Astrocyte'
    plot_regiontype_upset <- function(st, tr, ct){
        # create upset plot
        head(formatted_df)

        gene_lists <- formatted_df %>%
            filter(adj_pval < 0.05) %>%
            filter(formatted_st == st,
                   trait == tr,
                   formatted_ct == ct) %>%
            #filter(region_type == 'Full genes') %>%
            distinct(gene_id, formatted_st, formatted_rt) %>%
            unique() %>%
            group_by(formatted_st, formatted_rt) %>%
            summarise(genes=list(gene_id),
                      num_genes=n()
            ) %>%
            rename(summary_type=formatted_st, 
                   region_type=formatted_rt
            )
        head(gene_lists)
        gene_lists

        freqs <- gene_lists %>%
            arrange(-num_genes)
        head(freqs)
        ordered_genes <- as.character(freqs$region_type)
        ordered_genes


        # create upset plot for the same data
        selectlists <- gene_lists$genes
        names(selectlists) <- as.character(gene_lists$region_type) 
        names(selectlists)
        selectlists

        region_type_colors
        sub_colors <- region_type_colors[ordered_genes]
        length(sub_colors)
        sub_colors



        savefile <- paste0(sub_dir, tr, '_', st, "_", ct,  "_allRegions_celltype_upset.pdf")
        savefile
        pdf(savefile, width=7, height=7)
        p <- UpSetR::upset(fromList(selectlists),
              order.by='freq',
              nsets=10,
              keep.order=TRUE,
              nintersects=10,
              sets.bar.color=sub_colors,
              mb.ratio=c(0.6,0.4),
              text.scale=c(2, 2, 2, 2, 2, 2)
        )
        print(p)
        dev.off()
    }
    #plot_regiontype_upset(st='PCs', 'ceradsc', "Astro")
    plot_regiontype_upset(st='Avgs', 'ceradsc', "Astrocyte")



    get_counts <- function(){
        # filter fo significant p-vals
        # count the number of sig reults for each group
        head(formatted_df)
        sig_counts <- formatted_df %>%
            filter(adj_pval < 0.05) %>%
            select(gene_id,trait, summary_type, region_type, cell_type, full_array) %>%
            unique() %>%
            count(region_type, trait, summary_type, cell_type, full_array, name='sig_genes')
        head(sig_counts)


        # count the total number of input features for each group
        head(formatted_df)
        total_counts <- formatted_df %>%
            select(gene_id, trait, summary_type, region_type, cell_type, full_array) %>%
            unique() %>%
            count(region_type, trait, summary_type, cell_type, full_array, name='total_count') 
        head(total_counts)


        # compute the percentage of sig results for each group
        sig_counts <- sig_counts %>%
            full_join(total_counts) %>%
            mutate(percent = sig_genes / total_count * 100) %>%
            mutate(sig_genes = ifelse(is.na(sig_genes), 0, sig_genes)) %>%
            mutate(label = paste0(sig_genes)) 
        head(sig_counts)


        # change the region type for cpgs
        # set the position for text label
        # set color for summary types
        formatted_counts <- sig_counts %>%
            #filter(summary_type %in% c('PCs', 'Avgs')) %>%
            unique() %>%
            ungroup() %>%
            mutate(max_dmp = max(sig_genes)) %>%
            mutate(ypos = ifelse(sig_genes < max_dmp-(max_dmp/2),
                                 sig_genes,
                                 sig_genes-(max_dmp*(1/2)))) %>%
            mutate(color = ifelse(sig_genes < max_dmp-(max_dmp/2),
                                  summary_type, 'white')) 
        head(formatted_counts)

        list(sig_counts=sig_counts,
             formatted_counts=formatted_counts
        )
    }
    counts <- get_counts()
    sig_counts <- counts$sig_counts
    formatted_counts <- counts$formatted_counts

    sig_counts %>% head()

    compare_proportions <- function(){
        head(formatted_df)


        # compare proportions 
        head(formatted_counts)
        wide_counts <- formatted_counts %>%
            select(region_type, trait, summary_type, cell_type,
                   sig_genes, total_count) 
            #spread(summary_type, sig_genes) 
            #mutate(Avgs = as.numeric(Avgs)) %>%
            #mutate(PCs = as.numeric(PCs)) %>%
            #mutate(total_count = as.numeric(total_count)) 
        head(wide_counts)

        row <- wide_counts[1,]
        row

        ct = 'Oligo/OPC'
        st = 'Avgs'
        rt1 = '1to5kb'
        rt2 = 'Promoter'
        test_regiontype_props <- function(rt1, rt2, st, ct){
            print(paste("processing", rt1, rt2, st, ct))
            head(formatted_counts)
            test_counts <- formatted_counts %>%
                filter(summary_type == st,
                       cell_type == ct,
                       region_type %in% c(rt1, rt2)
                ) 
            test_counts
            rownames(test_counts) <- test_counts$region_type

            n1 <- test_counts[rt1,]$sig_genes
            n2 <- test_counts[rt2,]$sig_genes

            total1 <- test_counts[rt1,]$total_count
            total2 <- test_counts[rt2,]$total_count

            print(paste(rt1, n1, total1, round(n1/total1,4)))
            print(paste(rt2, n2, total2, round(n2/total2,4)))

            res <- prop.test(x=c(n1, n2),
                             n=c(total1, total2),
                      alternative='two.sided'
            )
            res
            names(res)
            res$conf.int[[1]]

            data.frame(region_type1 = rt1,
                       num_sig1 = n1,
                       total1 = total1,
                       region_type2 = rt2,
                       num_sig2 = n2,
                       total2 = total2,
                       pval = res$p.value,
                       lower_CI = res$conf.int[[1]],
                       upper_CI = res$conf.int[[2]]
            )
        }


        rts <- unique(formatted_counts$region_type) %>% as.character()
        rts <- setdiff(rts, 'Astro enh')
        rt2 <- 'Full gene'
        res <- lapply(rts, test_regiontype_props, rt2, st, ct)
        res_df <- do.call(rbind, res) %>%
            mutate(prop1 = num_sig1 / total1,
                   prop2 = num_sig2 / total2,
                   sig = ifelse(pval < 0.05, '*', ''),
                   greater = prop2 > prop1
            )
        res_df


        props <- apply(wide_counts, 1, test_props)
        props_df <- do.call(rbind, props)
        head(props_df)
        test_props <- function(row){
            #print(row)
            pcs_n <- row[['PCs']] %>% as.numeric()
            avgs_n <- row[['CpGs']] %>% as.numeric()
            total_n <- row[['total_count']] %>% as.numeric()
            #print(paste(pcs_n, avgs_n, total_n))

            res <- prop.test(x=c(pcs_n, avgs_n),
                             n=c(total_n, total_n),
                      alternative='greater'
            )

            data.frame(region_type=row[['region_type']],
                       trait=row[['trait']],
                       cell_type=row[['cell_type']],
                       pval=res$p.value)
        }
        props <- apply(wide_counts, 1, test_props)
        props_df <- do.call(rbind, props)
        head(props_df)

        props_corrected <- props_df %>%
            group_by(trait, greater_st, lesser_st, cell_type) %>%
            mutate(adj_pval = p.adjust(pval, method='bonferroni'))
        head(props_corrected)

        table(props_corrected$pval < 0.05)
        table(props_corrected$adj_pval < 0.05)
        head(wide_counts)

    }


    # plot the counts of significant DMPs
    with_cpgs=FALSE
    tr = "ceradsc"
    plot_full_gene_counts <- function(with_cpgs, tr){
        head(formatted_counts)

        filtered_counts <- formatted_counts %>%
            filter(summary_type != 'CpGs') %>%
            filter(region_type == 'full_gene', 
                   trait == tr
            ) 
        if (with_cpgs){
            filtered_counts <- formatted_counts %>%
                filter(region_type == 'Full gene', 
                       trait == tr
                ) 
        }

        head(filtered_counts)
        p <- filtered_counts %>%
            left_join(summary_type_map) %>%
            left_join(cell_type_map) %>%
            mutate(formatted_trait = "CERAD Score") %>%
            ggplot(aes(x=formatted_ct, y=sig_genes, fill=formatted_st)) +
            geom_bar(stat='identity', position='dodge') +
            geom_text(aes(label=label,
                          #y=ypos,
                          y=sig_genes,
                          #color=color
                          ),
                      hjust=0.5,
                      vjust=-0.5,
                      show.legend=FALSE,
                      size=5,
                      position=position_dodge(width=0.9)
                      ) +
            #facet_grid(region_type ~ cell_type, scales='fixed') +
            scale_color_manual(values=summary_type_colors, name="Summary type") + #, breaks=summary_types) +
            scale_fill_manual(values=summary_type_colors, name="Summary type") + #, breaks=summary_types) +
            theme_bw() +
            theme(
                  strip.background=element_rect(fill='white', color='black'),
                  panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  text=element_text(size=20, family='Open Sans'),
                  legend.title=element_text(size=20, face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
            ) +
            ggtitle(paste(tr, correction_type, "sig genes")) +
            xlab("Summary types") +
            ylab("Significant DM genes")

            main_filename = "_fullGene_dm_genes_counts_barplot.png"
            # not including cpgs
            savefile <- paste0(sub_dir, tr, "_", correction_type, 
                               main_filename)
            if (full_array){
                savefile <- paste0(sub_dir, tr, "_", correction_type, 
                                   "_fullArray", 
                                   main_filename)
            }
            if (with_cpgs){
                # inlucde cpgs
                savefile <- paste0(sub_dir, tr, "_", correction_type, 
                                   "_withCpgs",
                                   main_filename)
                if (full_array){
                    savefile <- paste0(sub_dir, tr, "_", correction_type, 
                                       "_withCpgs",
                                       "_fullArray", 
                                       main_filename)
                }
            }
            savefile

            ggsave(savefile, plot=p, width=12, height=7, units='in')


            main_filename = "_fullGene_dm_genes_counts_barplot.rds"
            savefile <- paste0(sub_dir, tr, "_", correction_type, 
                               main_filename)
            savefile
            saveRDS(p, savefile)
        1
    }

    #with_cpgs=FALSE
    #plot_full_gene_counts(with_cpgs=TRUE)
    #plot_full_gene_counts(with_cpgs=FALSE, 'ceradsc')


    # plot the counts of significant DMPs with traits
    with_cpgs=FALSE
    tr = 'ceradsc'
    rt = 'Full gene'
    plot_full_gene_counts_with_traits <- function(with_cpgs, rt, tr=""){
        head(formatted_counts)
        unique(formatted_counts$trait)

        filtered_counts <- formatted_counts %>%
            filter(summary_type != 'CpGs') %>%
            filter(region_type == rt) 
        if (with_cpgs){
            filtered_counts <- formatted_counts %>%
                filter(region_type == rt) 
        }

        if (!all(is.na(tr))){
            filtered_counts <- filtered_counts %>%
                filter(trait == tr)
        }

        # update the y position 
        head(filtered_counts)
        filtered_counts <- filtered_counts %>%
            group_by(trait) %>%
            mutate(max_dmp = max(sig_genes)) %>%
            mutate(ypos = sig_genes) %>%
            mutate(ypos = ifelse(sig_genes < max_dmp-(max_dmp/2),
                                 sig_genes,
                                 sig_genes-(max_dmp*(1/2)))) %>%
            mutate(color = ifelse(sig_genes < max_dmp-(max_dmp/2),
                                  as.character(summary_type), 'white'))  %>%
            #mutate(color = as.character(summary_type)) %>%
            mutate(formatted_trait = case_when(trait == 'cat_braaksc' ~ "Braak",
                                       trait == 'ceradsc' ~ "CERAD",
                                       trait == 'new_apoe' ~ "APOE",
                                       trait == "new_diag" ~ "Diagnosis"
                                       ) 
            )
        head(filtered_counts)

        summary_type_colors

        p <- filtered_counts %>%
            ggplot(aes(x=cell_type, y=sig_genes, fill=summary_type)) +
            geom_bar(stat='identity', position='dodge',
                     show.legend=TRUE) +
            geom_text(aes(label=label,
                          y=ypos,
                          color=color
                          ),
                      hjust=0.5,
                      vjust=-0.5,
                      show.legend=FALSE,
                      size=5,
                      position=position_dodge(width = .9)
                      ) +
            #facet_grid(formatted_trait ~ ., scales='free_y') +
            scale_color_manual(values=summary_type_colors, guide="none") +
            scale_fill_manual(name="Summary type", values=summary_type_colors) +
            #coord_cartesian(ylim=c(0,max(filtered_counts$ypos)+100)) +
            theme_bw() +
            theme(
                  strip.background=element_rect(fill='white', color='black'),
                  panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  text=element_text(size=20, family='Open Sans'),
                  legend.title=element_text(size=20, face='bold'),
                  legend.position='top',
                  legend.justification='right',
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
                  axis.text.x=element_text(angle=30, vjust=1, hjust=1)
            ) +
            #coord_cartesian(ylim=c(0,270)) +
            #ggtitle(paste(correction_type, "sig genes")) +
            xlab("Cell types") +
            #ylab(paste0(tr, "-associated genes"))
            ylab("Number of DM genes")

            main_filename = paste(tr, str_replace(tolower(rt), " ", "_"), "traits_fullGene_dm_genes_counts_barplot.png", sep='_')
            # not including cpgs
            savefile <- paste0(sub_dir,
                               main_filename)
            if (full_array){
                savefile <- paste0(sub_dir, correction_type, 
                                   "_fullArray", 
                                   main_filename)
            }
            if (with_cpgs){
                # inlucde cpgs
                savefile <- paste0(sub_dir, correction_type, 
                                   "_withCpgs",
                                   main_filename)
                if (full_array){
                    savefile <- paste0(sub_dir, correction_type, 
                                       "_withCpgs",
                                       "_fullArray", 
                                       main_filename)
                }
            }
            savefile

            #ggsave(savefile, plot=p, width=8, height=7, units='in')
            ggsave(savefile, plot=p, width=8, height=8, units='in')
        1
    }
    
    unique(formatted_df$region_type)
    #plot_full_gene_counts_with_traits(with_cpgs=FALSE, tr="")
    plot_full_gene_counts_with_traits(with_cpgs=TRUE, rt='Full gene', tr='ceradsc')
    plot_full_gene_counts_with_traits(with_cpgs=TRUE, rt='Promoter', tr='ceradsc')
    plot_full_gene_counts_with_traits(with_cpgs=TRUE, rt='Gene body', tr='ceradsc')


    # plot the counts of significant DMPs with traits
    with_cpgs=FALSE
    ct="Bulk"
    plot_allRegions_oneCelltype <- function(with_cpgs, ct, tr){
        head(formatted_counts)

        filtered_counts <- formatted_counts %>%
            filter(summary_type != 'CpGs') %>%
            filter(cell_type == ct)
        if (with_cpgs){
            filtered_counts <- formatted_counts %>%
                filter(cell_type == ct)
        }

        filtered_counts <- filtered_counts %>%
            filter(trait == tr)

        # update the y position 
        filtered_counts <- filtered_counts %>%
            #mutate(max_dmp = max(sig_genes)) %>%
            #mutate(ypos = ifelse(sig_genes < max_dmp-(max_dmp/2),
                                 #sig_genes,
                                 #sig_genes-(max_dmp*(1/2)))) %>%
            #mutate(color = ifelse(sig_genes < max_dmp-(max_dmp/2),
                                  #as.character(summary_type), 'white')) %>%
            mutate(ypos=sig_genes,
                  color=as.character(summary_type)
            )
        head(filtered_counts)

        summary_type_colors

        p <- filtered_counts %>%
            ggplot(aes(x=region_type, y=sig_genes, fill=summary_type)) +
            geom_bar(stat='identity', position='dodge',
                     show.legend=TRUE) +
            geom_text(aes(label=label,
                          y=ypos,
                          color=color
                          ),
                      hjust=0.5,
                      vjust=-0.5,
                      show.legend=FALSE,
                      size=5,
                      position=position_dodge(width = .9)
                      ) +
            facet_grid(trait ~ ., scales='fixed') +
            scale_color_manual(values=summary_type_colors, guide="none") +
            scale_fill_manual(name="Summary type", values=summary_type_colors) +
            theme_bw() +
            theme(
                  strip.background=element_rect(fill='white', color='black'),
                  panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  text=element_text(size=20, family='Open Sans'),
                  legend.title=element_text(size=20, face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
                  axis.text.x=element_text(angle=30, vjust=1, hjust=1)
            ) +
            #coord_cartesian(ylim=c(0,270)) +
            ggtitle(paste(correction_type, "sig genes")) +
            xlab("Cell types") +
            ylab("Significant DM genes")

            main_filename = paste0("_traits_allRegions_", ct, "_dm_genes_counts_barplot.png")
            # not including cpgs
            savefile <- paste0(sub_dir, correction_type, 
                               main_filename)
            if (full_array){
                savefile <- paste0(sub_dir, correction_type, 
                                   "_fullArray", 
                                   main_filename)
            }
            if (with_cpgs){
                # inlucde cpgs
                savefile <- paste0(sub_dir, correction_type, 
                                   "_withCpgs",
                                   main_filename)
                if (full_array){
                    savefile <- paste0(sub_dir, correction_type, 
                                       "_withCpgs",
                                       "_fullArray", 
                                       main_filename)
                }
            }
            savefile

            ggsave(savefile, plot=p, width=10, height=7, units='in')
        1
    }
    
    plot_allRegions_oneCelltype(with_cpgs=FALSE, ct="Bulk", tr="ceradsc")
    plot_allRegions_oneCelltype(with_cpgs=TRUE, ct="Bulk", tr="ceradsc")


    # plot with other region types
    with_cpgs=FALSE
    trait = 'allTraits'
    plot_region_gene_counts <- function(with_cpgs, trait){
        head(formatted_counts)

        t = trait
        filtered_counts <- formatted_counts %>%
             filter(trait == t) 
         if (t == ''){
            filtered_counts <- formatted_counts 
         }
        filtered_counts <- filtered_counts %>%
            left_join(summary_type_map) %>%
            left_join(region_type_map) %>%
            filter(summary_type != 'CpGs') %>%
            unique() %>%
            ungroup() %>%
            mutate(max_dmp = max(sig_genes)) %>%
            mutate(ypos = ifelse(sig_genes < max_dmp-(max_dmp/3),
                                 sig_genes,
                                 sig_genes-(max_dmp*(1/3)))) %>%
            mutate(color = ifelse(sig_genes < max_dmp-(max_dmp/3),
                                  formatted_st, 'white'))  
        if (with_cpgs){
            filtered_counts <- formatted_counts 
        }
        head(filtered_counts)

        trait_map <- data.frame(trait=unique(filtered_counts$trait)) 
        trait_map['formatted_tr'] <- c('Braak Stage', 'CERAD Score', 'APOE', 'Diagnosis')
        head(trait_map)

        p <- filtered_counts %>%
            left_join(summary_type_map) %>%
            left_join(region_type_map) %>%
            left_join(cell_type_map) %>%
            left_join(trait_map) %>%
            ggplot(aes(x=formatted_rt, y=sig_genes, fill=formatted_st)) +
            geom_bar(stat='identity', position='dodge') +
            geom_text(aes(label=label,
                          y=ypos,
                          color=color
                          ),
                      hjust=0.5,
                      vjust=-0.5,
                      show.legend=FALSE,
                      size=5,
                      position=position_dodge(width=0.9)
                      ) +
            facet_grid(formatted_ct ~ formatted_tr, scales='fixed') +
            scale_color_manual(values=summary_type_colors) + #, breaks=unique(filtered_counts$summary_type)) +
            scale_fill_manual(values=summary_type_colors, name="Summary type") +#, breaks=unique(filtered_counts$summary_types)) +
            theme_bw() +
            theme(
                  strip.background=element_rect(fill='white', color='black'),
                  panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  axis.text.x=element_text(angle=30, vjust=1, hjust=1),
                  text=element_text(size=20, family='Open Sans'),
                  legend.title=element_text(size=16, face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
            ) +
            #coord_cartesian(ylim=c(0,150)) +
            ggtitle(paste(trait, correction_type, "sig genes")) +
            xlab("Region types") +
            ylab("Number of significant DM genes")

            main_filename = "_allRegions_dm_genes_counts_barplot.png"
            # not including cpgs
            savefile <- paste0(sub_dir, trait, "_", correction_type, 
                               main_filename)
            if (full_array){
                savefile <- paste0(sub_dir, trait, "_", correction_type, 
                                   "_fullArray", 
                                   main_filename)
            }
            if (with_cpgs){
                # inlucde cpgs
                savefile <- paste0(sub_dir, trait, "_", correction_type, 
                                   "_withCpgs",
                                   main_filename)
                if (full_array){
                    savefile <- paste0(sub_dir, trait, "_", correction_type, 
                                       "_withCpgs",
                                       "_fullArray", 
                                       main_filename)
                }
            }
            savefile

        ggsave(savefile, plot=p, width=15, height=9, units='in')
        1
    }
    #plot_region_gene_counts(with_cpgs=TRUE)
    plot_region_gene_counts(with_cpgs=FALSE, trait='braaksc')
    plot_region_gene_counts(with_cpgs=FALSE, trait='ceradsc')
    #plot_region_gene_counts(with_cpgs=TRUE, trait='ceradsc')
    plot_region_gene_counts(with_cpgs=FALSE, trait='cat_braaksc')
    plot_region_gene_counts(with_cpgs=FALSE, trait='new_diag')
    plot_region_gene_counts(with_cpgs=FALSE, trait='new_apoe')


    # compare region types
    compare_region_types <- function(){
        # compare full gene results to the other region types within it
        head(formatted_df)
        sig_df <- formatted_df %>%
            filter(adj_pval < 0.05) %>%
            select(gene_id,trait, summary_type, region_type, cell_type, full_array) 
        head(sig_df)

        # get the sub regions of full gene
        unique(sig_df$region_type)
        contained_regions <- c("Promoters", "Exons", "Introns", "1to5kb", "5UTRs", "3UTRs")

        contained_df <- sig_df %>%
            filter(region_type %in% contained_regions) %>%
            select(-full_array, -region_type) %>%
            mutate(in_contained = TRUE)
        head(contained_df)

        full_df <- sig_df %>%
            filter(region_type == "Full genes") %>%
            select(-region_type, -full_array) %>%
            mutate(in_full = TRUE)
        head(full_df)

        st <- 'Avgs'
        combined_df <- full_df %>%
            full_join(contained_df) %>%
            filter(summary_type == st) %>%
            unique() %>%
            mutate(in_contained = ifelse(is.na(in_contained), FALSE, in_contained)) %>%
            mutate(in_full = ifelse(is.na(in_full), FALSE, in_full)) %>%
            mutate(dm_gene = ifelse(in_full & in_contained, "both", NA),
                   dm_gene = ifelse(!in_full & in_contained, "sub_only", dm_gene),
                   dm_gene = ifelse(in_full & !in_contained, "full_only", dm_gene)
            )
        head(combined_df)

        combined_df %>%
            filter(dm_gene == 'sub_only') %>%
            head()


        p <- ggplot(combined_df) +
            geom_bar(aes(x=cell_type, fill=dm_gene),
                     position='dodge'
            ) +
            facet_grid(trait ~ .) +
            theme_bw() +
            theme(
                  axis.text.x=element_text(angle=30, vjust=1, hjust=1),
                  text=element_text(size=20, family='Open Sans'),
                  #legend.title=element_text(size=16, face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
            ) +
            #coord_cartesian(ylim=c(0,150)) +
            ggtitle(paste(trait, st, correction_type, "sig genes")) +
            xlab("Cell types") +
            ylab("Number of significant DM genes")

        main_filename = "_compare_regions.png"
        # not including cpgs
        savefile <- paste0(sub_dir, st, correction_type, 
                           main_filename)
        if (full_array){
            savefile <- paste0(sub_dir, st, correction_type, 
                               "_fullArray", 
                               main_filename)
        }
        savefile

        ggsave(savefile, plot=p, width=8, height=7, units='in')
    }

    with_cpgs=FALSE
    plot_region_gene_counts_with_traits <- function(with_cpgs){
        head(formatted_counts)

        filtered_counts <- formatted_counts %>%
            filter(!is.na(cell_type)) %>%
            #filter(summary_type != 'CpGs') %>%
            unique() %>%
            ungroup() %>%
            mutate(max_dmp = max(sig_genes)) %>%
            mutate(ypos = ifelse(sig_genes < max_dmp-(max_dmp/2),
                                 sig_genes,
                                 sig_genes-(max_dmp*(1/2)))) %>%
            mutate(color = ifelse(sig_genes < max_dmp-(max_dmp/2),
                                  summary_type, 'white')) 
        if (with_cpgs){
            filtered_counts <- formatted_counts 
        }

        p <- filtered_counts %>%
            ggplot(aes(x=cell_type, y=sig_genes, fill=summary_type)) +
            geom_bar(stat='identity', position='dodge') +
            geom_text(aes(label=label,
                          y=ypos,
                          color=color
                          ),
                      hjust=0.5,
                      vjust=-0.5,
                      show.legend=FALSE,
                      size=5,
                      ) +
            facet_grid(trait ~ region_type, scales='fixed') +
            scale_color_manual(values=summary_type_colors, breaks=summary_types) +
            scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
            theme_bw() +
            theme(
                  axis.text.x=element_text(angle=80, vjust=0.5),
                  text=element_text(size=16, family='Open Sans'),
                  legend.title=element_text(size=16, face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
            ) +
            #coord_cartesian(ylim=c(0,150)) +
            ggtitle(paste(trait, correction_type, "sig genes")) +
            xlab("Summary types") +
            ylab("Number of significant DM genes")

            main_filename = "_allRegions_withTraits_dm_genes_counts_barplot.png"
            # not including cpgs
            savefile <- paste0(sub_dir, trait, "_", correction_type, 
                               main_filename)
            if (full_array){
                savefile <- paste0(sub_dir, trait, "_", correction_type, 
                                   "_fullArray", 
                                   main_filename)
            }
            if (with_cpgs){
                # inlucde cpgs
                savefile <- paste0(sub_dir, trait, "_", correction_type, 
                                   "_withCpgs",
                                   main_filename)
                if (full_array){
                    savefile <- paste0(sub_dir, trait, "_", correction_type, 
                                       "_withCpgs",
                                       "_fullArray", 
                                       main_filename)
                }
            }
            savefile

        ggsave(savefile, plot=p, width=15, height=9, units='in')
        1
    }
    #plot_region_gene_counts_with_traits(with_cpgs=FALSE)

    # plot the counts of significant DMPs using a subset of region types
    restrict_region_types <- function(){
        head(formatted_counts)
        p <- formatted_counts %>%
            filter(region_type %in% c('Full genes', 'Promoters')) %>%
            #filter(cell_type %in% c('Bulk')) %>%
            ggplot(aes(x=summary_type, y=sig_genes, fill=summary_type)) +
            geom_bar(stat='identity', position='dodge') +
            geom_text(aes(label=label,
                          y=ypos,
                          color=color
                          ),
                      hjust=0.5,
                      vjust=-0.5,
                      show.legend=FALSE,
                      size=11,
                      ) +
            facet_grid(region_type ~ cell_type, scales='fixed') +
            scale_color_manual(values=summary_type_colors, breaks=summary_types) +
            scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
            theme_bw() +
            theme(
                  strip.background=element_rect(fill='white', color='black'),
                  panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  text=element_text(size=20, family='Open Sans'),
                  legend.title=element_text(size=20, face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
            ) +
            ggtitle(paste(trait, correction_type, "sig genes")) +
            xlab("Summary types") +
            ylab("Number of significant DM genes")+
            coord_cartesian(ylim=c(0,100))


        savefile <- paste0(sub_dir, trait, "_", correction_type, 
                           "_fullgene_dm_genes_counts_barplot.png")
        if (full_array){
            savefile <- paste0(sub_dir, trait, "_", correction_type, 
                               "_fullgene_full_array_dm_genes_counts_barplot.png")
        }
        savefile

        ggsave(savefile, plot=p, width=18, height=8, units='in')
        ggsave(savefile, plot=p, width=4, height=8, units='in')
        1
    }


    plot_percentages <- function(){
        # plot the percentages of the counts
        head(formatted_counts)
        ymax = max(formatted_counts$percent)
        p <- formatted_counts %>%
            filter(full_array = TRUE) %>%
            mutate(label=paste0(round(percent, 2), "%")) %>%
            ggplot(aes(x=summary_type, y=percent, fill=summary_type)) +
            geom_bar(stat='identity', position='dodge') +
            geom_text(aes(label=label,
                          y=percent,
                          ),
                      color='black',
                      hjust=0.5,
                      vjust=-0.5,
                      show.legend=FALSE,
                      size=8,
                      position=position_dodge(width=1)
                      ) +
            facet_grid(region_type ~ cell_type, scales='fixed') +
            #scale_color_manual(values=summary_type_colors, breaks=summary_types) +
            scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
            theme_bw() +
            theme(
                  text=element_text(size=20, family='Open Sans'),
                  legend.title=element_text(size=20, face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
            ) +
            ggtitle(paste(trait, correction_type, "Sig DMPs")) +
            xlab("Summary types") +
            ylab("Percentage of significant DM regions") +
            coord_cartesian(ylim=c(0,ymax+0.2))

        savefile <- paste0(sub_dir, trait, "_", correction_type, 
                           "_dm_genes_percentage_barplot.png")
        if (full_array){
            savefile <- paste0(sub_dir, trait, "_", correction_type, 
                               "_full_array_dm_genes_percentage_barplot.png")
        }
        savefile

        ggsave(savefile, plot=p, width=18, height=13, units='in')
        1
    }
    #plot_percentages()



    plot_total_features <- function(){
        ## -- plot the number of total features per summary type
        head(sig_counts)
        counts <- sig_counts %>%
            group_by(region_type) %>%
            rename(total_probes=total_count) %>%
            mutate(max_dmp=max(total_probes)) %>%
            mutate(ypos=ifelse(total_probes < max_dmp-(max_dmp/2),
                               total_probes,
                               total_probes-max_dmp*(1/2))) %>%
            mutate(color=ifelse(total_probes < max_dmp-(max_dmp/2),
                                as.character(summary_type),
                                'black')) %>%
            select(-cell_type) %>%
            unique()
        head(counts)

        summary_type_colors['black'] <- 'black'

        p <- counts %>%
            ggplot(aes(x=summary_type, y=total_probes, fill=summary_type)) +
            geom_bar(stat='identity', position='dodge') +
            geom_text(aes(label=total_probes,
                          y=ypos,
                          color=color
                          ),
                      hjust=0.5,
                      vjust=-0.5,
                      show.legend=FALSE,
                      size=10
                      #position=position_dodge(width=1)
                      ) +
            facet_grid(region_type ~ ., scales='free_y') +
            scale_color_manual(values=summary_type_colors, breaks=summary_types) +
            scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
            theme_bw() +
            theme(
                  text=element_text(size=20, family='Open Sans'),
                  legend.title=element_text(size=20, face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0))
            ) +
            ggtitle(paste(trait, correction_type, "number of genes")) +
            xlab("Summary types") +
            ylab("Number of total genes")

        savefile <- paste0(sub_dir, trait, "_", correction_type, 
                           "_gene_count_barplot.png")
        if (full_array){
            savefile <- paste0(sub_dir, trait, "_", correction_type, 
                               "_full_array_gene_count_barplot.png")
        }
        savefile

        ggsave(savefile, plot=p, width=18, height=10, units='in')
    }


    ## -- compare the genes to the opentargets genes
    load_ad_genes <- function(){
        # load the opentargets genes
        datafile <- paste0("/home/eulalio/deconvolution/rosmap/data/gene_sets/AD_opentargets.csv")
        opentargets <- read_csv(datafile)
        head(opentargets)

        # map the gene_ids in our dmp results to symbols
        #datafile <- paste0("/home/eulalio/deconvolution/data/gene_biomart_table.csv")
        #gene_annots <- read_csv(datafile) %>%
            #rename(gene_no_vers=gene_id,
                   #gene_id=ensembl_gene_id_version) 
        #head(gene_annots)

        head(formatted_df)
        head(opentargets)
        #head(gene_annots)

        # add gene annotations to the opentargets genes
        opentargets_annots <- opentargets %>%
            select(symbol, overallAssociationScore)
        head(opentargets_annots)

        # join the opentargets genes with the dm genes
        head(formatted_df)
        genes_novers <- unique(formatted_df$gene_id) %>%
            str_remove("\\.[0-9]+")
        head(genes_novers)

        savefile <- paste0(savedir, "gene_symbols.rds")
        if (!file.exists(savefile)){
            ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
            dm_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                     #filters = "ensembl_gene_id", 
                     #values = unique(genes_novers), 
                     mart = ensembl) 
            saveRDS(dm_symbols, savefile)
        } else{
            print("Loading existing biomart genes file")
            dm_symbols <- readRDS(savefile)
        }
        head(dm_symbols)

        dm_symbols_formatted <- dm_symbols %>%
            rename(gene_novers=ensembl_gene_id,
                   symbol=hgnc_symbol
            )

        length(unique(genes_novers))
        dim(dm_symbols)

        setdiff(genes_novers, dm_symbols$ensembl_gene_id)


        head(formatted_df)
        head(opentargets_annots)
        dm_opentargets <- formatted_df %>%
            mutate(gene_novers = str_remove(gene_id, "\\.[0-9]+")) %>%
            left_join(dm_symbols_formatted) %>%
            left_join(opentargets) %>%
            mutate(sig_dm=adj_pval < 0.05)
        head(dm_opentargets)

        dim(dm_opentargets)

        # load the NIAGDS genes also
        datafile <- paste0("/home/eulalio/deconvolution/data/gene_sets/NIAGDS_AD_genes_list.csv")
        niagds <- read_csv(datafile, col_names=c('symbol'))
        head(niagds)

        head(dm_opentargets)
        dm_opentargets <- dm_opentargets %>%
            mutate(niagds = symbol %in% niagds$symbol)



        list(dm_opentargets=dm_opentargets,
             opentargets=opentargets_annots)
    }
    res <- load_ad_genes()
    all_dm_opentargets=res$dm_opentargets
    opentargets=res$opentargets
    head(all_dm_opentargets)


    dm_opentargets=all_dm_opentargets %>%
        filter(sig_dm)
    head(dm_opentargets)

    visualize_groups <- function(){
        dm_counts <- all_dm_opentargets %>%
            filter(summary_type == 'PCs',
                   adj_pval < 0.05) %>%
            count(gene_id, summary_type, region_type, cell_type, trait) %>%
            filter(n > 1)
        head(dm_counts)

        tmp <- all_dm_opentargets %>%
            inner_join(dm_counts) %>%
            arrange(-overallAssociationScore)
        head(tmp)
        dim(tmp)
        unique(tmp$symbol)
    }

    # create glm contrast data
    create_glm_data <- function(){
        head(all_dm_opentargets)

        sub_data <- all_dm_opentargets %>%
            select(gene_id, symbol, overallAssociationScore, adj_pval, summary_type, region_type, cell_type) %>%
            mutate(ad_gene = as.numeric(!is.na(overallAssociationScore))) %>%
            mutate(summary_type = paste(summary_type, "DM", sep='_')) %>%
            group_by(summary_type, region_type, cell_type, gene_id) %>%
            arrange(adj_pval) %>%
            slice(1) %>%
            mutate(DM = as.numeric(adj_pval < 0.05)) %>%
            select(-adj_pval) %>%
            unique() %>%
            spread(summary_type, DM, fill=0)
        head(sub_data)

        savefile <- paste0(savedir, "bh_all_ad_dm_glm_formatted_results.rds")
        savefile
        saveRDS(sub_data, savefile)
        

        filtered_data <- sub_data %>%
            filter(region_type == 'Full gene',
                   cell_type == 'Astrocyte'
            ) %>%
            select(gene_id, ad_gene, Avgs_DM, PCs_DM) %>%
            unique() 
            #remove_rownames() %>%
            #column_to_rownames('gene_id')
        head(filtered_data)

        savefile <- paste0(savedir, "bh_full_gene_astro_ad_dm_glm_formatted_results.rds")
        savefile
        saveRDS(filtered_data, savefile)

        table(filtered_data$ad_gene)
        table(filtered_data$Avgs_DM)
        table(filtered_data$PCs_DM)


        
    }


    quant=0.95
    tr = 'ceradsc'
    plot_focus_genes <- function(quant, tr){
        ## -- focus on just a subset of the genes and region types
        head(dm_opentargets)
        #summary(opentargets$overallAssociationScore)
        #quant=0.7
        (thresh=quantile(opentargets$overallAssociationScore, quant))

        head(all_dm_opentargets)
        c <- 0.000001
        sub_dm <- all_dm_opentargets %>%
            #filter(!is.na(overallAssociationScore)) %>%
            filter(formatted_rt == 'Full gene',
                   trait == 'ceradsc') %>%
            #filter(cell_type == 'Bulk') %>%
            mutate(full_score = ifelse(is.na(overallAssociationScore), 0, overallAssociationScore)) %>%
            mutate(minus_p = 1-P.Value) %>%
            mutate(logit_pval = log((minus_p+c) / (1-minus_p+c)),
                   scaled_pval = scale(logit_pval)) %>%
            mutate(logit_score = log((full_score+c) / (1-full_score+c)),
                   scaled_score=scale(logit_score)
            ) 
        head(sub_dm)
        summary(sub_dm$full_score)

        unique(sub_dm$region_type)
        save_data <- sub_dm %>%
            select(summary_type, cell_type, gene_id, logit_pval, logit_score) %>%
            filter(summary_type == 'avgs') %>%
            select(-summary_type)
        head(save_data)

        (savefile <- paste0(savedir, "avgs_ceradsc_fullGene_dmLogitPvals_adScores.csv"))
        write_csv(save_data, savefile)        



        #res <- lm(full_score ~ bh_pval + cell_type + summary_type, data=sub_dm)
        #res <- lm(logit_score ~ logit_pval + cell_type + summary_type, data=sub_dm)
        #res <- lm(scaled_score ~ scaled_pval*summary_type + cell_type, data=sub_dm)
        res <- lm(scaled_score ~ cell_type + scaled_pval*summary_type + 0, data=sub_dm)
        res <- lm(logit_score ~ cell_type + logit_pval*summary_type + 0, data=sub_dm)

        res <- lm(logit_score ~ logit_pval*summary_type + cell_type+0, data=sub_dm)
        summary(res)

        # checking lrt 
        st_dm <- sub_dm %>%
            filter(summary_type == 'avgs') %>%
            filter(cell_type == 'oligo_opc')
            #mutate(cell_type = factor(cell_type, levels=c('bulk', 'astro', 'endo', 'neuron', 'oligo_opc')))
        head(st_dm)
        unique(st_dm$cell_type)
        summary(lm(logit_score ~ logit_pval, data=st_dm))
        res1 <- lm(logit_score ~ cell_type * logit_pval, data=st_dm)
        res2 <- lm(logit_score ~ cell_type + logit_pval, data=st_dm)
        summary(res1)
        summary(res2)

        anova(res1, res2, test='LRT')

        #logit_test <- data.frame(x=seq(0,1,0.01)) %>%
            #mutate(logit_x=log((x+c) / (1-x+c)))
        #head(logit_test)

        # up arrow
        up_arrow_position <- data.frame(x=min(sub_dm$logit_pval), y=-9.5)
        up_text_position <- data.frame(x=-10.5, y=-9.5-0.1)

        # side arrow
        side_arrow_position <- data.frame(x=14, y=-11)
        side_text_position <- data.frame(x=4, y=-10.85)

        head(sub_dm)
        min_score=min(sub_dm$overallAssociationScore, na.rm=TRUE)
        (savefile <- paste0(sub_dir, trait, "_full_gene_pval_adScore_scatter.png"))
        p <- sub_dm %>%
            #filter(!is.na(overallAssociationScore)) %>%
            ggplot(aes(y=logit_score, x=logit_pval)) +
            geom_smooth(method='lm', se=FALSE, linewidth=2, aes(color=formatted_ct, linetype=formatted_st)) +
            #geom_point(aes(color=formatted_st, shape=formatted_st)) +
            # draw up arrow
            geom_segment(data=up_arrow_position, aes(x=x, y=y-0.25, xend=x, yend=y),
                         arrow=arrow(length=unit(0.3, "cm")), 
                         linetype='solid',
                         color='black',
                         linewidth=1.1
            ) +
            geom_text(data=up_text_position, aes(x=x, y=y, label="More\nevidence"), 
                      color='black',
                      hjust=0,
                      size=6
                      ) +
            # draw side arrow
            geom_segment(data=side_arrow_position, aes(x=x-5, y=y, xend=x, yend=y),
                         arrow=arrow(length=unit(0.3, "cm")), 
                         linetype='solid',
                         color='black',
                         linewidth=1.1
            ) +
            geom_text(data=side_text_position, aes(x=x, y=y, label="Higher\nsignificance"), 
                      color='black',
                      hjust=0,
                      size=6
                      ) +
            theme_bw() +
            theme(
                  strip.background=element_rect(fill='white', color='black'),
                  #panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  text=element_text(size=18),
                  legend.key.width = unit(2, "line")
                  #legend.position='top',
                  #legend.box='vertical'
                  #axis.text.x=element_text(angle=0, hjust=0, vjust=1, size=12)
            ) +
            scale_color_manual(values=cell_type_colors, name="Cell type") +
            scale_linetype_discrete(name='Summary type')+
            ylab("logit(AD Association Score)") +
            xlab("logit(1 - DM p-value)")
        savefile
        ggsave(p, file=savefile, width=7, height=5)

        (savefile <- paste0(sub_dir, trait, "_full_gene_pval_adScore_scatter.rds"))
        saveRDS(p, savefile)


        head(dm_opentargets)
        filtered_genes <- dm_opentargets %>%
            filter(region_type == 'full_gene') %>%
            select(gene_id, trait, symbol, 
                   summary_type=formatted_st, region_type=formatted_rt, cell_type=formatted_ct, 
                   overallAssociationScore, niagds) %>%
            unique() %>%
            filter(summary_type != 'CpGs') %>%
            filter(overallAssociationScore > thresh | niagds) %>%
            filter(trait == tr) 
        head(filtered_genes)

        length(unique(filtered_genes$symbol))


        # handles case where one summary type didn't find anything
        # filter thse out later
        filtered_genes2 <- rbind(filtered_genes, c('TEST', NA, NA, 'Avgs', NA, NA, 0, FALSE))
        filtered_genes2 <- rbind(filtered_genes2, c('TEST', NA, NA, 'rPCs', NA, NA, 0, FALSE))

        filtered_genes %>%
            filter(symbol == "PSENEN")


        head(filtered_genes2)
        wide_genes <- filtered_genes2 %>%
            mutate(hit = TRUE) %>%
            spread(summary_type, hit, fill=FALSE) %>%
            mutate(hit = ifelse(Avgs & !rPCs, 'Avgs only', 'rPCs only'),
                   hit = ifelse(rPCs & Avgs, 'Both', hit)
            ) %>%
            filter(gene_id != 'TEST') %>%
            mutate(overallAssociationScore = as.numeric(overallAssociationScore))
        head(wide_genes)
        tail(wide_genes)
        
        prop.test(x=c(1, 26), n=c(348, 348), alternative='less')

        wide_genes %>%
            filter(!is.na(overallAssociationScore)) %>%
            #filter(Avgs, !rPCs) %>%
            filter(!Avgs, rPCs) %>%
            select(symbol) %>%
            unique() %>%
            nrow()

        focus_genes <- wide_genes %>%
            filter(!is.na(overallAssociationScore)) %>%
            mutate(symbol = ifelse(niagds, paste0("*", symbol), symbol))
        head(focus_genes)
        tail(focus_genes)

        pc_only <- focus_genes %>%
            filter(!Avgs, PCs) 
            #arrange(-overallAssociationScore)
        head(pc_only)
        length(unique(pc_only$gene_id))
        

        #summary_colors <- c(brewer.pal(n=3, name='Dark2')[1:2], 'grey')
        #summary_colors
        #summary_type_colors
        summary_colors <- c(summary_type_colors[c('Avgs', 'rPCs')], '#99e6e6')
        summary_colors

        hit_colors <- setNames(summary_colors, c('Avgs only', 'rPCs only', 'Both'))


        # plot tileplot for paper
        head(focus_genes)
        p <- focus_genes %>%
            filter(region_type %in% c('Full gene')) %>%
            select(gene_id, symbol, cell_type, overallAssociationScore, hit) %>%
            spread(cell_type, hit) %>%
            gather('cell_type', 'hit', Bulk:"Oligo/OPC") %>%
            ggplot() +
            geom_tile(aes(x=cell_type,
                          y=reorder(symbol, overallAssociationScore),
                          fill=hit),
                      color='black',
                      linewidth=1
                      ) +
            #facet_grid(trait ~ region_type, scales='free_y') +
            coord_flip() +
            theme_bw() +
            theme(
                  #axis.text.y=element_blank(),
                  axis.ticks.x=element_blank(),
                  strip.background=element_rect(fill='white', color='black'),
                  panel.spacing=unit(0, 'lines'),
                  #panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  panel.border=element_blank(),
                  text=element_text(size=30, family='Open Sans'),
                  axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                  axis.title.y=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank()
            ) +
            scale_fill_manual(values=hit_colors, breaks=c('Avgs only', 'rPCs only', 'Both'),
                              name="Differential\nmethylation",
                              na.value='white'
                              ) +
            ylab(paste("Genes ordered by AD Association Score")) 
            #xlab(paste("Cell type")) +
            #ggtitle(paste("AD Open Targets significant hits, quantile", quant,
                          #"pval correction", correction_type, tr))

        savefile <- paste0(sub_dir, tr, "_", correction_type, 
                           "_quant", round(quant,2),
                           "_opentargets_focusgenes_tileplot.png")
        if (full_array){
            savefile <- paste0(sub_dir, tr, "_", correction_type, 
                               "_quant", round(quant,2),
                               "_full_array_opentargets_focusgenes_tileplot.png")
        }
        savefile

        ggsave(savefile, plot=p, width=17, height=5, units='in')

        savefile <- paste0(sub_dir, tr, "_", correction_type, 
                           "_quant", round(quant,2),
                           "_opentargets_focusgenes_tileplot.rds")
        saveRDS(p, savefile)

        # plot tileplot for slides
        head(focus_genes)
        p <- focus_genes %>%
            filter(region_type %in% c('Full gene')) %>%
            select(gene_id, symbol, cell_type, overallAssociationScore, hit, niagds) %>%
            spread(cell_type, hit) %>%
            gather('cell_type', 'hit', Bulk:"Oligo/OPC") %>%
            mutate(cell_type = fct_relevel(cell_type, 'Bulk', 'Astrocyte')) %>%
            ggplot() +
            geom_tile(aes(x=cell_type,
                          y=reorder(symbol, overallAssociationScore),
                          fill=hit),
                      color='black',
                      linewidth=1
                      ) +
            facet_wrap(. ~ niagds, scales='free_y') +
            #facet_grid(trait ~ region_type, scales='free_y') +
            #coord_flip() +
            theme_bw() +
            theme(
                  #axis.text.y=element_blank(),
                  axis.ticks.x=element_blank(),
                  #strip.background=element_rect(fill='white', color='white'),
                  strip.text.x=element_blank(),
                  panel.spacing=unit(0, 'lines'),
                  #panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  panel.border=element_blank(),
                  text=element_text(size=30, family='Open Sans'),
                  axis.text.x=element_text(angle=30, hjust=1, vjust=1),
                  axis.title.x=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank()
            ) +
            scale_fill_manual(values=hit_colors, breaks=c('Avgs only', 'rPCs only', 'Both'),
                              name="Differential\nmethylation",
                              na.value='white'
                              ) +
            ylab(paste("Genes ordered by\nAD Association Score")) 
            #xlab(paste("Cell type")) +
            #ggtitle(paste("AD Open Targets significant hits, quantile", quant,
                          #"pval correction", correction_type, tr))

        savefile <- paste0(sub_dir, tr, "_", correction_type, 
                           "_quant", round(quant,2),
                           "_opentargets_focusgenes_tileplot_for_slides.png")
        savefile
        ggsave(savefile, plot=p, width=14, height=8, units='in')


        # plot all region types
        head(focus_genes)
        top_n <- focus_genes %>%
            select(symbol, overallAssociationScore) %>%
            unique() %>%
            arrange(-overallAssociationScore) %>%
            top_n(15)
        top_n

        p <- focus_genes %>%
            filter(symbol %in% top_n$symbol) %>%
            ggplot() +
            geom_tile(aes(x=cell_type,
                          y=reorder(symbol, overallAssociationScore),
                          fill=hit)) +
            facet_grid(trait ~ region_type, scales='free_y') +
            theme_bw() +
            theme(
                  #axis.text.y=element_blank(),
                  #axis.ticks.y=element_blank(),
                  strip.background=element_rect(fill='white', color='black'),
                  panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  text=element_text(size=20, family='Open Sans'),
                  axis.text.x=element_text(angle=60, hjust=1, vjust=1)
            ) +
            scale_fill_manual(values=hit_colors, breaks=c('Avgs only', 'PCs only', 'Both'),
                              name="Differential\nmethylation") +
            ylab(paste("Genes ordered by AD Association Score")) +
            xlab(paste("Cell type")) +
            ggtitle(paste("AD Open Targets significant hits, quantile", quant,
                          "pval correction", correction_type, tr))

        savefile <- paste0(sub_dir, tr, "_", correction_type, 
                           "_allRegions",
                           "_quant", round(quant,2),
                           "_opentargets_focusgenes_tileplot.png")
        if (full_array){
            savefile <- paste0(sub_dir, tr, "_", correction_type, 
                               "_quant", round(quant,2),
                               "_full_array_opentargets_focusgenes_tileplot.png")
        }
        savefile

        ggsave(savefile, plot=p, width=16, height=9, units='in')



        # plot counts of AD genes found
        head(filtered_genes)
        tail(filtered_genes)
        counts <- filtered_genes %>%
            filter(!is.na(overallAssociationScore)) %>%
            count(summary_type, region_type, cell_type) %>%
            spread(cell_type, n, fill=0) %>%
            gather('cell_type', 'n', -summary_type, -region_type) %>%
            group_by(region_type) %>%
            mutate(max=max(n)) %>%
            mutate(ypos=ifelse(n < max-(max/2),
                               n,
                               n-max*(1/2))) %>%
            mutate(color=ifelse(n < max-(max/2),
                                as.character(summary_type),
                                'black')) 
        head(counts)


        p <- counts %>%
            #filter(region_type %in% c('Full genes', 'Promoters')) %>%
            ggplot(aes(x=region_type, y=n, fill=factor(summary_type, levels=c('PCs', 'Avgs')))) +
            geom_bar(stat='identity', position='dodge') +
            #geom_text(aes(label=n,
                          #y=ypos,
                          #),
                      #hjust=0.5,
                      #vjust=-0.5,
                      #show.legend=FALSE,
                      #size=10,
                        #color='black'
                      ##position=position_dodge(width=1)
                      #) +
            facet_grid(cell_type ~ ., scales='fixed') +
            scale_fill_manual(values=summary_type_colors, breaks=c('PCs', 'Avgs'),
                              name="Summary\ntype", drop=FALSE) +
            #scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
            scale_x_discrete(drop=FALSE) +
            theme_bw() +
            theme(
                  strip.background=element_rect(fill='white', color='black'),
                  panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  text=element_text(size=20, family='Open Sans'),
                  legend.title=element_text(size=20),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
                  axis.text.x=element_text(angle=30, hjust=1, vjust=0.9)
            ) +
            ggtitle(paste(tr, correction_type, "number of genes")) +
            xlab("Cell types") +
            ylab("Number of AD genes")

        savefile <- paste0(sub_dir, tr, "_", correction_type, 
                           "_quant", round(quant,2),
                           "_opentargets_focusgenes_barplot.png")
        if (full_array){
            savefile <- paste0(sub_dir, tr, "_", correction_type, 
                               "_quant", round(quant,2),
                               "_full_array_opentargets_focusgenes_barplot.png")
        }
        savefile

        ggsave(savefile, plot=p, width=13, height=8, units='in')
    }
    quant=0.95
    plot_focus_genes(quant=quant, 'ceradsc')
    plot_focus_genes(quant=quant, 'cat_braaksc')
    plot_focus_genes(quant=quant, 'new_diag')
    plot_focus_genes(quant=quant, 'new_apoe')
    #plot_focus_genes(quant=quant, 'braaksc')


    quant=0
    compare_ad_proportions <- function(quant){
        # compare proportions 
        head(opentargets)
        (thresh=quantile(opentargets$overallAssociationScore, quant))

        head(all_dm_opentargets)

        (savefile <- paste0(sub_dir, trait, "_DM_AD_counts.png"))
        p <- all_dm_opentargets %>%
            mutate(ad_gene = !is.na(overallAssociationScore)) %>%
            mutate(group = case_when(sig_dm & ad_gene ~ "DM, AD",
                                     !sig_dm & ad_gene ~ "-DM, AD",
                                     sig_dm & !ad_gene ~ "DM, -AD",
                                     !sig_dm & !ad_gene ~ "-DM, -AD"
                                     )) %>%
            filter(group %in% c('DM, AD', 'DM, -AD')) %>%
            ggplot() +
            geom_bar(aes(x=group, fill=summary_type),
                     position=position_dodge()) +
            facet_grid(cell_type ~ region_type, scales='free_y') +
            scale_fill_manual(values=summary_type_colors, name="Summary type") +
            theme_bw() +
            theme(
                  strip.background=element_rect(fill='white', color='black'),
                  #panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
                  text=element_text(size=14, family='Open Sans'),
                  legend.title=element_text(face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
                  axis.text.x=element_text(angle=50, vjust=1, hjust=1)
            ) 
        ggsave(p, file=savefile, width=11, height=8)



        ad_counts <- all_dm_opentargets
        runs <- expand.grid(cell_type = unique(ad_counts$cell_type),
                       region_type = unique(ad_counts$region_type),
                       trait = unique(ad_counts$trait),
                       st1 = unique(ad_counts$summary_type),
                       st2 = unique(ad_counts$summary_type)
        ) %>%
            filter(st1 != st2) %>%
            mutate(rn = row_number())
        head(runs)
        runs %>%
            filter(region_type == 'Full gene')

        row <- runs[49,]
        compare_st_props <- function(row){
            ct <- row[['cell_type']] 
            rt <- row[['region_type']] %>% as.character()
            tr <- row[['trait']] %>% as.character()
            st1 <- row[['st1']] %>% as.character()
            st2 <- row[['st2']] %>% as.character()
            rn <- row[['rn']] %>% as.character()

            print(paste(ct, rt, tr, st1, st2, rn))

            # filter down to data for this run
            head(ad_counts)
            run_counts <- ad_counts %>%
                filter(trait == tr, 
                       cell_type == ct,
                       region_type == rt
                ) %>% 
                # make sure these are DM genes only
                filter(adj_pval < 0.05) %>%
                select(gene_novers, symbol, overallAssociationScore, summary_type) %>%
                #unique() %>%
                # add AD or not column
                mutate(thresholded_score = ifelse(overallAssociationScore < thresh, NA, overallAssociationScore)) %>%
                mutate(ad_gene = ifelse(is.na(thresholded_score), 'not_AD', 'AD'))
            head(run_counts)

            if (nrow(run_counts) == 0){
                return(NA)
            }
            if (length(unique(run_counts$summary_type)) < 2){
                return(NA)
            }

            # create a contingency table
            (tab <- table(run_counts$ad_gene, run_counts$summary_type))
            tab[1,1] / sum(tab[,1])
            tab[1,2] / sum(tab[,2])

            # run fisher's exact test
            (res <- chisq.test(tab))
            (res <- fisher.test(tab))

            names(res)
            res$alternative
            tab

            tab[1,st1]

            data.frame(region_type = rt,
                       cell_type = ct,
                       summary_type1 = st1,
                       summary_type2 = st2,
                       pval = res$p.value,
                       lower_CI = res$conf.int[[1]],
                       upper_CI = res$conf.int[[2]],
                       odds_ratio = res$estimate[[1]],
                       st1_ad = tab[['AD', st1]],
                       st2_ad = tab[['AD', st2]],
                       st1_nad = tab[['not_AD', st1]],
                       st2_nad = tab[['not_AD', st2]]
            ) 
        }

        res <- apply(runs, 1, compare_st_props)
        res_df <- do.call(rbind, res)
        head(res_df)



        formatted_res <- res_df %>%
            mutate(sig = pval < 0.05) %>%
            mutate(prop1 = st1_ad / (st1_ad + st1_nad)) %>%
            mutate(prop2 = st2_ad / (st2_ad + st2_nad)) 
        head(formatted_res)

        formatted_res %>%
            filter(region_type == 'Full gene')

        savefile <- paste0(savedir, "summary_type_dm_ad_proportions_fishers_test_results.rds")
        savefile
        saveRDS(formatted_res, savefile)
        

        unique(all_dm_opentargets$region_type)

        head(formatted_df)
        head(dm_opentargets)

        ad_total_counts <- all_dm_opentargets %>%
            select(gene_novers, overallAssociationScore, summary_type, cell_type, region_type, trait) %>%
            mutate(overallAssociationScore = ifelse(overallAssociationScore < thresh, NA, overallAssociationScore)) %>%
            mutate(ad_gene = !is.na(overallAssociationScore)) %>%
            unique() %>%
            count(summary_type, region_type, cell_type, trait, ad_gene) %>%
            group_by(summary_type, region_type, cell_type, trait) %>%
            mutate(total_n = sum(n)) %>%
            filter(ad_gene)
        head(ad_total_counts)

        # count total and ad genes
        head(dm_opentargets)
        ad_dm_counts <- all_dm_opentargets %>%
            filter(adj_pval < 0.05) %>%
            mutate(overallAssociationScore = ifelse(overallAssociationScore < thresh, NA, overallAssociationScore)) %>%
            mutate(ad_gene = !is.na(overallAssociationScore)) %>%
            select(gene_id, summary_type, region_type, cell_type, trait, ad_gene) %>%
            unique() %>%
            count(summary_type, region_type, cell_type, trait, ad_gene) %>%
            group_by(summary_type, region_type, cell_type, trait) %>%
            mutate(total_n = sum(n)) %>%
            filter(ad_gene)
        head(ad_dm_counts)

        head(ad_total_counts)
        ad_total_props <- ad_total_counts %>%
            rename(total_ad=n, total_genes=total_n) %>%
            mutate(total_prop=total_ad/total_genes)
        head(ad_total_props)

        ad_dm_props <- ad_dm_counts %>%
            rename(sig_ad=n, total_sig=total_n) %>%
            mutate(sig_prop=sig_ad/total_sig)
        head(ad_dm_props)

        joined_props <- ad_total_props %>%
            left_join(ad_dm_props) %>%
            arrange(region_type, cell_type)
        head(joined_props)

        joined_props[is.na(joined_props)] <- 0
        joined_props %>%
            filter(region_type == '5UTR')

        ad_counts <- joined_props
        runs <- expand.grid(cell_type = unique(ad_counts$cell_type),
                       region_type = unique(ad_counts$region_type),
                       trait = unique(ad_counts$trait),
                       st = unique(ad_counts$summary_type)
                       #st2 = unique(ad_counts$summary_type)
        ) #%>%
            #filter(st1 != st2)
        head(runs)
        dim(runs)
        runs

        row <- runs[1,]
        compare_total_props <- function(row){
            ct <- row[['cell_type']] 
            rt <- row[['region_type']] %>% as.character()
            tr <- row[['trait']] %>% as.character()
            #st1 <- row[['st1']] %>% as.character()
            #st2 <- row[['st2']] %>% as.character()
            st <- row[['st']] %>% as.character()

            print(paste(ct, rt, tr, st))

            #head(ad_counts)
            wide_counts <- ad_counts %>%
                filter(cell_type == ct,
                       region_type == rt,
                       trait == tr,
                       summary_type == st
                ) 

            if (nrow(wide_counts) == 0){
                return(NA)
            }

            head(wide_counts)
            n1_label = paste0("total_ad")
            n2_label = paste0("sig_ad")
            t1_label = paste0("total_genes")
            t2_label = paste0("total_sig")

            n1 <- wide_counts[[1,n1_label]] %>% as.numeric()
            n2 <- wide_counts[[1,n2_label]] %>% as.numeric()
            t1 <- wide_counts[[1,t1_label]] %>% as.numeric()
            t2 <- wide_counts[[1,t2_label]] %>% as.numeric()

            print(paste("total", n1, t1, round(wide_counts$total_prop,2), "sig", n2, t2, round(wide_counts$sig_prop,2)))

            if (wide_counts$sig_prop == 0){
                print("no sig ad genes")
                return(NA)
            }

            #res <- tryCatch(
                #expr={
                    #prop.test(x=c(n1, n2),
                         #n=c(t1, t2),
                          #alternative='greater')},
                #error=function(cond){
                    #print('err') 
                    #return(NA)},
                #warning=function(cond){
                    #print('warn') 
                    #return(NA)}
            #)
            res <- prop.test(x=c(n1, n2),
                             n=c(t1, t2),
                      alternative='less'
            )
            res
            names(res)
            res$conf.int[[1]]

            data.frame(region_type = rt,
                       cell_type = ct,
                       summary_type = st,
                       ad_total = n1,
                       genes_total = t1,
                       ad_sig = n2,
                       genes_sig = t2,
                       pval = res$p.value,
                       lower_CI = res$conf.int[[1]],
                       upper_CI = res$conf.int[[2]]
            ) %>%
            mutate(
                       prop_total = ad_total / genes_total,
                       prop_sig = ad_sig / genes_sig,
                       greater = prop_sig > prop_total
                       )
        }
        props <- apply(runs , 1, compare_total_props)
        props_df <- do.call(rbind, props) %>%
            na.omit()
        head(props_df)
        dim(props_df)

        table(props_df$pval < 0.05)

        props_corrected <- props_df %>%
            group_by(region_type, summary_type) %>%
            mutate(adj_pval = p.adjust(pval, method='BH')) %>%
            mutate(sig = ifelse(adj_pval < 0.05, '*', ''))
        head(props_corrected)

        props_corrected %>%
            filter(cell_type == 'Neuron',
                   region_type == 'Shelf'
            )

        table(props_corrected$pval < 0.05)
        table(props_corrected$adj_pval < 0.05)

        head(props_corrected)
        long_props <- props_corrected %>%
            select(region_type:summary_type, prop_total, prop_sig) %>%
            spread(summary_type, prop_sig, fill=0) %>%
            rename(Total=prop_total) %>%
            gather('summary_type', 'prop', Total, Avgs:PCs) %>%
            left_join(props_corrected[c('region_type', 'cell_type', 'summary_type', 'adj_pval')]) %>%
            mutate(adj_pval = ifelse(summary_type == 'Total', 1, adj_pval)) %>%
            mutate(sig = ifelse(adj_pval < 0.1, '.', ''),
                   sig = ifelse(adj_pval < 0.05, '*', sig),
                   sig = ifelse(adj_pval < 0.001, '**',sig)
            ) %>%
            unique() 
        head(long_props)

        foldchange <- long_props %>%
            select(-adj_pval, -sig) %>%
            spread(summary_type, prop) %>%
            gather('summary_type', 'prop', Avgs:PCs) %>%
            mutate(fold_change = prop / Total) %>%
            left_join(long_props) %>%
            mutate(bin_sig = ifelse(sig == "", FALSE, TRUE))
        head(foldchange)

        long_props %>%
            filter(
                   region_type == 'Shelf',
                   cell_type == 'Neuron'
            )

        tmp_st_colors <- summary_type_colors
        i=8
        tmp_st_colors['Total'] = brewer.pal(i, 'Dark2')[i]
        head(foldchange)
        p <- ggplot(foldchange, aes(x=cell_type)) +
            geom_point(aes(y=fold_change, color=summary_type, size=bin_sig)) +
            geom_hline(yintercept=1, linetype="dotted") +
            #geom_bar(aes(fill=summary_type),
                     #stat='identity',
                     #position='dodge'
            #) +
            #geom_text(aes(label=sig, group=summary_type, y=prop+0.05, color=summary_type), 
                      #size=5,
                      #position = position_dodge(width = 0.9),
                      #) +
           facet_wrap(. ~ region_type, nrow=2) +
           scale_color_manual(values=tmp_st_colors, name="Summary type") +
           scale_fill_manual(values=tmp_st_colors, name="Summary type") +
           coord_flip() +
           theme_bw() +
           theme(text = element_text(size=20),
                 axis.text.x = element_text(angle=60, hjust=1, vjust=1)
           ) +
            xlab("Cell type") +
            ylab("Proportion of AD genes") +
            ggtitle(paste("Proportion of significant AD genes", quant, 'ceradsc'))


        (savefile <- paste0(sub_dir, tr, "_", correction_type, 
                           "_quant", round(quant,2),
                           "_opentargets_proportion_dotplot.png"))
        ggsave(savefile, plot=p, width=13, height=8, units='in')


        # can we run a fisher's test against background 
        head(dm_opentargets)

        runs <- expand.grid(cell_type = unique(dm_opentargets$cell_type),
                       region_type = unique(dm_opentargets$region_type),
                       trait = unique(dm_opentargets$trait),
                       st = unique(dm_opentargets$summary_type)
        ) %>%
        mutate(rn = row_number())
        head(runs)
        dim(runs)
        runs

        row <- runs[1,]

        sub_data <- all_dm_opentargets %>%
            select(gene_id, symbol, overallAssociationScore, adj_pval, summary_type, region_type, cell_type) %>%
            mutate(ad_gene = as.numeric(!is.na(overallAssociationScore))) %>%
            mutate(summary_type = paste(summary_type, "DM", sep='_')) %>%
            group_by(summary_type, region_type, cell_type, gene_id) %>%
            arrange(adj_pval) %>%
            slice(1) %>%
            mutate(DM = as.numeric(adj_pval < 0.05)) %>%
            select(-adj_pval) %>%
            unique() %>%
            spread(summary_type, DM, fill=0)
        head(sub_data)

        savefile <- paste0(savedir, "bh_all_ad_dm_glm_formatted_results.rds")
        savefile
        saveRDS(sub_data, savefile)

        sub_data <- readRDS(savefile)
        head(sub_data)


        formatted_dat <- sub_data %>%
            gather('summary_type', 'DM', Avgs_DM:PCs_DM) %>%
            mutate(summary_type = str_remove(summary_type, "_DM")) %>%
            mutate(pc = as.numeric(summary_type == 'PCs')) %>%
            mutate(full_score = ifelse(is.na(overallAssociationScore), 0, overallAssociationScore))
        head(formatted_dat)

        pc_genes <- formatted_dat %>%
            filter(summary_type == 'PCs')
        head(pc_genes)

        avg_genes <- formatted_dat %>%
            filter(summary_type == 'Avgs')


        full_gene_dat <- formatted_dat %>%
            filter(region_type == 'Full gene')
        head(full_gene_dat)

        full_lm <- lm(full_score ~ cell_type + DM*summary_type, data=full_gene_dat)
        summary(full_lm)





        # plot the proportions of ad for dm and not dm genes across
        # summary types, cell types, region types
        counts <- formatted_dat %>%
            group_by(region_type, cell_type, summary_type, DM) %>%
            summarize(ad_prop = sum(ad_gene) / sum(ad_gene %in% c(0,1))) %>%
            mutate(DM = as.boolean(DM))
            #spread(DM, ad_prop) 
        head(counts)


        head(all_dm_opentargets)
        total_ad_prop <- all_dm_opentargets %>%
            select(gene_id, overallAssociationScore) %>%
            unique() %>%
            mutate(ad_gene = !is.na(overallAssociationScore)) %>%
            summarise(ad_prop = sum(ad_gene) / (sum(ad_gene) + sum(!ad_gene)))
        total_ad_prop

        p <- ggplot(counts) +
            geom_bar(aes(x=summary_type, y=ad_prop, fill=as.factor(DM)),
                     stat='identity', position=position_dodge(width=0.9)
            ) +
            facet_grid(region_type ~ cell_type)
        (savefile <- paste0(sub_dir, "full_ad_props.png"))
        ggsave(p, file=savefile)

        # create line plot
        head(counts)
        p <- counts %>%
            filter(region_type == 'Full gene') %>%
            ggplot() +
            geom_hline(yintercept=total_ad_prop[[1]], linetype='dashed') +
            geom_point(aes(x=cell_type, y=ad_prop, color=summary_type, fill=summary_type, shape=as.factor(DM), group=summary_type),
                       alpha=0.9, position=position_dodge(width=0.6), size=3
                       ) +
            #facet_wrap(. ~ cell_type, ncol=1) +
            scale_color_manual(values=summary_type_colors, name="Summary type") +
            scale_fill_manual(values=summary_type_colors, name="Summary type") +
            #scale_shape_manual(values=c(1,2)) +
            #scale_size_manual(values=c(2,3)) +
            #coord_cartesian(ylim=c(0,0.9)) +
            guides(shape=guide_legend(title="Differential\nmethylation")) +
            theme_bw() +
            theme(
                  strip.background=element_rect(fill='white', color='black'),
                  #panel.spacing=unit(0, 'lines'),
                  panel.border=element_rect(color='black', fill=NA, linewidth=0.9),
                  text=element_text(size=20, family='Open Sans'),
                  legend.title=element_text(size=20, face='bold'),
                  plot.title=element_text(hjust=0.5),
                  axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
                  axis.text.x=element_text(angle=30, vjust=1, hjust=1)
            ) +
            ggtitle(paste("CERAD score AD proportion")) +
            ylab("Proportion of Alzheimer's disease genes") +
            xlab("Region type")
        (savefile <- paste0(sub_dir, "full_ad_props_line.png"))
        ggsave(p, file=savefile, width=8, height=7)

        set.seed(1174117171)
        formatted_dat$permuted_st <- sample(formatted_dat$summary_type, size=nrow(formatted_dat), replace=FALSE)
        head(formatted_dat)
        table(formatted_dat$summary_type, formatted_dat$permuted_st)

        glm_dat <- formatted_dat %>%
            filter(region_type == 'Full gene')

        fit_glmer <- function(glm_dat){
            # fit model
            # fit with no intercept
            res <- glmer(ad_gene ~  pc * DM + (1 | region_type) + (1 | cell_type), data=formatted_dat, family='binomial')
            #savefile <- paste0(savedir, "glmer_res.rds")

            # fit with intercept
            res <- glmer(ad_gene  ~ pc + summary_type:DM + (1 | region_type) + (1 | cell_type), data=formatted_dat, family='binomial')
            savefile <- paste0(savedir, "glmer_res_nointercept.rds")

            # nicely formatted output
            res <- glm(ad_gene ~  cell_type + summary_type + summary_type:DM + 0, data=glm_dat, family='binomial')
            summary(res)
            (tidy_res <- tidy(res, conf.int=TRUE, exponentiate=TRUE))
            tidy_res


            # use a simpler model with no interaction term
            full_gene_simple_dat <- sub_data %>%
                filter(region_type == 'Full gene')
            head(full_gene_simple_dat)

            res <- glm(ad_gene ~ cell_type + Avgs_DM + PCs_DM + 0, data=full_gene_simple_dat, family='binomial')
            summary(res)
            (tidy_res <- tidy(res, conf.int=TRUE, exponentiate=TRUE))

            tidy_res
            formatted_res <- tidy_res %>%
                mutate(term = str_remove(term, 'cell_type')) %>%
                mutate(term = ifelse(grepl('DM', term), paste0("DM ", str_remove(term, "_DM")), term)) %>%
                mutate(formatted_cat = ifelse(grepl('DM', term), 'Summary type', 'Cell type'))
            formatted_res

            # plot these results
            #(savefile <- paste0(sub_dir, trait, "_glm_AD_prop_odds_ratios.png"))
            (savefile <- paste0(sub_dir, trait, "_lm_AD_prop_odds_ratios.png"))
            p <- formatted_res %>%
                ggplot(aes(y=reorder(term, estimate), x=estimate, color=formatted_cat)) +
                geom_point(size=3) +
                geom_linerange(aes(xmin=conf.low, xmax=conf.high)) +
                geom_vline(xintercept=1, linetype='dashed', color='grey') +
                theme_bw() +
                guides(color=guide_legend(title='Term type')) +
                scale_color_manual(values=c('cornflowerblue', 'firebrick1')) +
                theme(
                      text=element_text(size=20, family='Open Sans'),
                      legend.title=element_text(size=20, face='bold'),
                      plot.title=element_text(hjust=0.5)
                      #axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
                      #axis.text.x=element_text(angle=30, vjust=1, hjust=1)
                ) +
                xlab("Odds ratio") +
                ylab("Model term")
            ggsave(p, file=savefile, height=6, width=8)


            saveRDS(res, savefile)
            res <- readRDS(savefile)

            summary(res)
            (tidy_res <- tidy(res, conf.int=TRUE, exponentiate=TRUE))

            ranef(res)


            # model with no intercept, but less nicely formatted output
            res_noint_noformat <- glmer(ad_gene  ~ summary_type*DM + (1 | region_type) + (1 | cell_type), data=formatted_dat, family='binomial')

            savefile <- paste0(savedir, "glmer_res_nointercept_noformat.rds")
            saveRDS(res_noint_noformat, savefile)

            summary(res_noint_noformat)
            tidy(res_noint_noformat)


            # try permuting summary type 
            res_perm <- glmer(ad_gene  ~  permuted_st + permuted_st:DM + (1 | region_type) + (1 | cell_type) + 0, data=formatted_dat, family='binomial')
            savefile <- paste0(savedir, "glmer_res_permuted_summaryType.rds")
            savefile
            saveRDS(res_perm, savefile)

            summary(res_perm)
            (tidy_res <- tidy(res_perm, conf.int=TRUE, exponentiate=TRUE))
        }


        head(runs)
        runs %>%
            filter(cell_type == 'Neuron', 
                   region_type == '5UTR'
            )
        row <- runs[23,]
        compare_total_props <- function(row){
            ct <- row[['cell_type']] 
            rt <- row[['region_type']] %>% as.character()
            tr <- row[['trait']] %>% as.character()
            st <- row[['st']] %>% as.character()
            rn <- row[['rn']] %>% as.character

            print(paste(ct, rt, tr, st, rn, "out of", nrow(runs)))

            # subset the data
            head(sub_data)
            run_data <- sub_data %>%
                filter(cell_type == ct,
                       region_type == rt
                ) 
            head(run_data)

            # pivot data
            long_data <- run_data %>%
                gather('summary_type', 'DM', Avgs_DM:PCs_DM) %>%
                mutate(summary_type=str_remove(summary_type, "_DM")) %>%
                # filter for summary type
                filter(summary_type == st)
            head(long_data)

            tab <- table(long_data$ad_gene, long_data$DM)
            tab

            if (ncol(tab) < 2) return(NA)
            if (nrow(tab) < 2) return(NA)

            result <- fisher.test(tab, conf.level=0.95)
            or <- (result$estimate)
            ci <- result$conf.int
            lci <-  ci[[1]]
            uci <-  ci[[2]]
            pval <- result$p.value

            res <- data.frame(summary_type=st,
                              region_type=rt,
                              cell_type=ct,
                              odds_ratio=or,
                              lower_ci=lci,
                              upper_ci=uci,
                              pval=pval
            )
        }

        res <- apply(runs, 1, compare_total_props)
        res_df <- do.call(rbind, res)
        head(res_df)

        formatted_df <- res_df %>%
            group_by(summary_type) %>%
            mutate(adj_pval = p.adjust(pval, method='BH')) %>%
            mutate(log_or = log(odds_ratio),
                   log_lower_ci = log(lower_ci),
                   log_upper_ci = log(upper_ci)
            ) %>%
            na.omit() %>%
            mutate(sig = adj_pval < 0.05)
            #mutate(sig = pval < 0.05)
        head(formatted_df)        

        table(formatted_df$pval < 0.05)
        table(formatted_df$adj_pval < 0.05)

        head(formatted_df)
        p <- formatted_df %>%
            ggplot(aes(x=cell_type, y=log_or, color=summary_type)) +
            geom_point(aes(shape=sig),
                       size=3,
                       position=position_dodge(width=0.5)) +
            geom_linerange(aes(ymin=log_lower_ci, ymax=log_upper_ci),
                           position=position_dodge(width=0.5)) +
            geom_hline(yintercept=0, linetype="dotted") +
           facet_wrap(. ~ region_type, nrow=2, scales='fixed') +
           scale_color_manual(values=tmp_st_colors, name="Summary type") +
           scale_fill_manual(values=tmp_st_colors, name="Summary type") +
           coord_flip(ylim=c(-2,2)) +
           theme_bw() +
           theme(text = element_text(size=20),
                 axis.text.x = element_text(angle=0, hjust=0, vjust=0)
           ) +
            xlab("Cell type") +
            ylab("Fisher's Test Log Odds Ratio") +
            ggtitle(paste("Log Odds Ratio of DM AD genes", quant, 'ceradsc'))

        (savefile <- paste0(sub_dir, tr, "_", correction_type, 
                           "_quant", round(quant,2),
                           "_opentargets_logodds_dotplot.png"))
        ggsave(savefile, plot=p, width=13, height=8, units='in')
    }

    1
}



check_pcs <- function(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=TRUE){
    head(combined_df)


    print(paste("correction type", correction_type))

    full_array_str <- "full_array"
    if (!full_array) full_array_str <- "matched_array"

    sub_dir <- paste0(plotting_dir, correction_type, "_", full_array_str, "/")
    dir.create(sub_dir, showWarnings=FALSE)

    # set the p-value correction type
    combined_df$adj_pval <- combined_df$bh_pval
    if (correction_type == 'bf') combined_df$adj_pval <- combined_df$bf_pval


    # filter for using the full array or not
    if (full_array){
        combined_df <- combined_df %>%
            filter(full_array)
    } else{
        combined_df <- combined_df %>%
            filter(!full_array)
    }

    # find genes that are DM with PCs but not Avgs
    head(combined_df)

    sub_df <- combined_df %>%
        filter(formatted_rt == 'Full gene') %>%
        select(feature, trait:adj_pval) %>%
        separate(feature, c('gene_id', 'pc'), sep='-', fill='right') %>%
        mutate(sig_dm = adj_pval < 0.05) %>%
        group_by(gene_id, formatted_st, trait, formatted_ct, formatted_rt) %>%
        summarise(min_pval = min(adj_pval)) %>%
        mutate(sig_dm =  min_pval < 0.05) %>%
        select(-min_pval) %>%
        spread(formatted_st, sig_dm) %>%
        filter(!Avgs, PCs)
    head(sub_df)
    dim(sub_df)

    # how many pcs for each gene?
    head(combined_df)
    pc_counts <- combined_df %>%
        filter(formatted_rt == 'Full gene') %>%
        filter(summary_type == 'pcs') %>%
        select(feature) %>%
        separate(feature, c('gene_id', 'pc'), sep='-', fill='right') %>%
        count(gene_id, name='num_pcs')
    head(pc_counts)


    # how many cpgs per gene?
    datafile <- paste0(datadir, "full_gene_astro_gene_cpg_map.rds")
    cpg_map <- readRDS(datafile)
    head(cpg_map)

    cpg_count <- cpg_map %>%
        count(gene_id, name='num_cpgs')
    head(cpg_count)

    # add pc counts to selected genes
    select_genes <- sub_df %>%
        left_join(pc_counts) %>%
        left_join(cpg_count)
    head(select_genes)

    select_genes %>%
        filter(num_cpgs == 17
        ) 

    head(combined_df)
    sorted_df <- combined_df %>%
        filter(formatted_rt == 'Full gene',
               summary_type == 'pcs') %>%
        separate(feature, c('gene_id', 'pc'), sep='-', fill='right') %>%
        filter(gene_id %in% select_genes$gene_id) %>%
        arrange(-abs(logFC))
    head(sorted_df)
    

    combined_df %>%
        filter(grepl("ENSG00000076662.10", feature))

    pcs_df <- combined_df %>%
        filter(summary_type == 'pcs')
    head(pcs_df)


    # count number of sig PCs 
    formatted_pcs <- pcs_df %>%
        filter(adj_pval < 0.05) %>%
        select(feature, summary_type, region_type, cell_type, trait) %>%
        separate(feature, c('gene_id', 'pc'), sep='-')
    head(formatted_pcs)

    table(formatted_pcs$pc)

    # how many PCs are sig per gene?
    total_counts <- formatted_pcs %>%
        count(region_type, cell_type, trait, gene_id, name='sig_pcs')
    head(total_counts)

    combined_dat <- formatted_pcs %>%
        left_join(total_counts)
    head(combined_dat)
    unique(combined_dat$region_type)

    head(combined_dat)
    p <- combined_dat %>%
        filter(region_type == 'Full genes',
               trait == 'ceradsc') %>%
        mutate(sig_pcs = factor(sig_pcs)) %>%
        mutate(pc = factor(pc, levels=paste0("PC", 1:16))) %>%
        select(-pc) %>%
        unique() %>%
        ggplot(aes(x=cell_type, fill=sig_pcs)) +
        geom_bar(position='stack') +
        facet_grid(trait ~ ., scales='free_y') +
        #scale_color_manual(values=summary_type_colors, breaks=summary_types) +
        #scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
        theme_bw() +
        theme(
              strip.background=element_rect(fill='white', color='black'),
              panel.spacing=unit(0, 'lines'),
              panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
              text=element_text(size=20, family='Open Sans'),
              legend.title=element_text(size=20, face='bold'),
              plot.title=element_text(hjust=0.5),
              axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
              axis.text.x=element_text(angle=30, vjust=1, hjust=1)
        ) +
        #scale_y_continuous(trans='log10') +
        ggtitle(paste(correction_type, "sig genes")) +
        ylab("Count") +
        xlab("Significant PCs per gene")
    savefile <- paste0(sub_dir, "combined_pc_counts_bar.png")
    savefile
    ggsave(p, file=savefile, width=7, height=7)


    # which pc is usually DM?
    head(formatted_pcs)
    pc_counts <- formatted_pcs %>%
        count(pc, summary_type, region_type, cell_type, trait, name='pc_count')
    head(pc_counts)


    head(combined_dat)
    ordered_pcs <- combined_dat %>%
        arrange(-sig_pcs) 
    (savefile <- paste0(sub_dir, "ceradsc_full_gene_dm_pcs_per_gene.csv"))
    write_csv(ordered_pcs, savefile)

    dim(combined_dat)
    p <- combined_dat %>%
        filter(region_type == 'Full gene',
               trait == 'ceradsc') %>%
        mutate(sig_pcs = factor(sig_pcs)) %>%
        mutate(pc = factor(pc, levels=paste0("PC", 1:16))) %>%
        ggplot(aes(x=cell_type, fill=pc)) +
        geom_bar(position='stack') +
        facet_wrap(. ~ sig_pcs, scales='free_y') +
        #scale_color_manual(values=summary_type_colors, breaks=summary_types) +
        #scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
        theme_bw() +
        theme(
              strip.background=element_rect(fill='white', color='black'),
              panel.spacing=unit(0, 'lines'),
              panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
              text=element_text(size=20, family='Open Sans'),
              legend.title=element_text(size=20, face='bold'),
              plot.title=element_text(hjust=0.5),
              axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
              axis.text.x=element_text(angle=50, vjust=1, hjust=1)
        ) +
        #scale_y_continuous(trans='log10') +
        ggtitle(paste(correction_type, "sig genes")) +
        ylab("Count") +
        xlab("Cell type")
    savefile <- paste0(sub_dir, "combined_pc_counts_bar2.png")
    savefile
    ggsave(p, file=savefile, width=9, height=7)

    # how many sig pcs per gene?
    counts <- formatted_pcs %>%
        count(region_type, cell_type, trait, gene_id, name='sig_pcs') %>%
        arrange(-sig_pcs) 
    head(counts)

    table(counts$sig_pcs)
    summary(counts$sig_pcs)


    combined_counts <- total_counts %>%
        select(-gene_id) %>%
        left_join(pc_counts) %>%
        unique()
    head(combined_counts)

    head(combined_dat)
    #p <- combined_counts %>%
    p <- combined_dat %>%
        filter(region_type == 'Full genes',
               trait == 'ceradsc') %>%
        #mutate(sig_pcs = as.factor(n)) %>%
        ggplot(aes(x=cell_type, fill=pc)) +
        #geom_bar(stat='identity', position='dodge') +
        #geom_text(aes(label=nn,
                      #y=nn,
                      ##color=color
                      #),
                  #hjust=0.5,
                  #vjust=-0.5,
                  #show.legend=FALSE,
                  #size=5,
                  ##position=position_dodge(width = .9)
                  #) +
        facet_grid(trait ~ ., scales='fixed') +
        #scale_color_manual(values=summary_type_colors, breaks=summary_types) +
        #scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
        theme_bw() +
        theme(
              strip.background=element_rect(fill='white', color='black'),
              panel.spacing=unit(0, 'lines'),
              panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
              text=element_text(size=20, family='Open Sans'),
              legend.title=element_text(size=20, face='bold'),
              plot.title=element_text(hjust=0.5),
              axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
        ) +
        #coord_cartesian(ylim=c(0,270)) +
        ggtitle(paste(correction_type, "sig genes")) +
        ylab("Count") +
        xlab("Significant PCs per gene")
    savefile <- paste0(sub_dir, "pcs_sig_per_gene_bar.png")
    savefile
    ggsave(p, file=savefile, width=12, height=7)

    # plto pc counts
    head(pc_counts)
    p <- pc_counts %>%
        filter(region_type == 'Full genes',
               trait == 'ceradsc'
        ) %>%
        mutate(pc = str_remove(pc, "PC")) %>%
        mutate(pc = factor(pc, levels=1:16)) %>%
        ggplot(aes(x=cell_type, y=pc_count, fill=pc),
               show.legend=FALSE) +
        geom_bar(stat='identity', position='stack') +
        #geom_text(aes(label=pc_count,
                      #y=pc_count,
                      ##color=color
                      #),
                  #hjust=0.5,
                  #vjust=-0.5,
                  #show.legend=FALSE,
                  #size=5,
                  ##position=position_dodge(width = .9)
                  #) +
        facet_grid(trait ~ ., scales='fixed') +
        #scale_color_manual(values=summary_type_colors, breaks=summary_types) +
        #scale_fill_manual(values=summary_type_colors, breaks=summary_types) +
        theme_bw() +
        theme(
              strip.background=element_rect(fill='white', color='black'),
              panel.spacing=unit(0, 'lines'),
              panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
              text=element_text(size=20, family='Open Sans'),
              legend.title=element_text(size=20, face='bold'),
              plot.title=element_text(hjust=0.5),
              axis.title.y=element_text(margin=margin(t=0, r=30, b=0, l=0)),
              axis.text.x=element_text(angle=30, vjust=1, hjust=1)
        ) +
        #coord_cartesian(ylim=c(0,270)) +
        ggtitle(paste(correction_type, "sig genes")) +
        ylab("Count") +
        xlab("Cell type")
    savefile <- paste0(sub_dir, "pcs_count.png")
    savefile
    ggsave(p, file=savefile, width=7, height=7)

}

plot_manhattan <- function(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=FALSE){
    head(combined_df)

    tr <- "ceradsc"
    st <- "PCs"
    rt <- "Full genes"
    ct <- "Astrocytes"

    create_manhattan_plot <- function(row){
        tr <- row[['trait']] %>% as.character()
        st <- row[['summary_type']] %>% as.character()
        rt <- row[['region_type']] %>% as.character() 
        ct <- row[['cell_type']] %>% as.character()

        filtered_df <- combined_df %>%
            filter(trait == tr,
                   summary_type == st,
                   region_type == rt,
                   cell_type == ct
            )

        # figure out x-axis positions
        # use TSS ofr the genes
        datafile <- paste0("/home/eulalio/deconvolution/data/region_annotations/hg38_genes_promoters_annotations.R")
        gene_annots <- read_tsv(datafile)
        head(gene_annots)

        # get fields to annotate genes
        formatted_annots <- gene_annots %>%
            select(chrom=seqnames, pos=end, gene_id=gencode_gene_id, symbol=gencode_gene_name)
        head(formatted_annots)

        head(filtered_df)
        dm_data <- filtered_df %>%
            select(feature, pval=bh_pval) %>%
            separate(feature, c("gene_id", "pc"), sep='-', fill='right') %>%
            left_join(formatted_annots) %>%
            mutate(chrom = as.numeric(str_remove(chrom, "chr"))) 
        head(dm_data)

        data_cum <- dm_data %>%
            group_by(chrom) %>%
            summarise(max_bp = max(pos)) %>%
            mutate(bp_add = lag(cumsum(max_bp), default=0)) %>%
            select(chrom, bp_add)
        head(data_cum)


        load_ad_genes <- function(){
            # load the opentargets genes
            datafile <- paste0("/home/eulalio/deconvolution/rosmap/data/gene_sets/AD_opentargets.csv")
            opentargets <- read_csv(datafile)
            head(opentargets)


            # map the gene_ids in our dmp results to symbols
            datafile <- paste0("/home/eulalio/deconvolution/data/gene_biomart_table.csv")
            gene_annots <- read_csv(datafile) %>%
                rename(gene_no_vers=gene_id,
                       gene_id=ensembl_gene_id_version) 
            head(gene_annots)

            head(opentargets)
            head(gene_annots)


            # add gene annotations to the opentargets genes
            opentargets_annots <- opentargets %>%
                left_join(gene_annots) %>%
                select(symbol, overallAssociationScore, gene_id)
            head(opentargets_annots)

            opentargets_annots
        }
        opentargets_annots <- load_ad_genes()
        head(opentargets_annots)

        (thresh <- quantile(opentargets_annots$overallAssociationScore, 0.8))
        opentargets_top <- opentargets_annots %>%
            filter(overallAssociationScore > thresh)

        head(dm_data)
        dm_annots <- dm_data %>%
            inner_join(data_cum, by="chrom") %>%
            mutate(bp_cum = pos + bp_add) %>%
            # add a highlight label
            mutate(is_highlight = ifelse(pval < 0.05, "yes", "no")) %>%
            mutate(is_annotate = ifelse((pval < 0.05) & (gene_id %in% opentargets_top$gene_id), "yes", "no")) %>%
            # add highly significant ones too
            mutate(is_annotate = ifelse(-log10(pval) > 3, "yes", is_annotate))
        head(dm_annots)

        dm_annots %>%
            filter(symbol == "PSENEN")

        table(dm_annots$is_highlight)
        table(dm_annots$is_annotate)

        # get the center of chromosome to position x axis labels
        axis_set <- dm_annots %>%
            group_by(chrom) %>%
            summarize(center = mean(bp_cum))
        head(axis_set)

        head(dm_annots)
        ylim <- -log10(min(dm_annots$pval, na.rm=TRUE)) + 2
        ylim

        sig <- 0.05

        manhplot <- ggplot(dm_annots,
                            aes(x=bp_cum, y=-log10(pval),
                                color=as_factor(chrom), size=-log10(pval))) +
            geom_hline(yintercept = -log10(sig), color="grey40", linetype="dashed") +
            geom_point(alpha = 0.75) +
            scale_x_continuous(label = axis_set$chrom, breaks=axis_set$center) +
            scale_y_continuous(expand = c(0,0), 
                               limits=c(0,ylim)
                               ) +
            scale_size_continuous(range = c(0.5,3)) +

            scale_color_manual(values = rep(c("darkred", "dodgerblue3"), unique(length(axis_set$chrom)))) +
            # add highlighted points
            geom_point(data=subset(dm_annots, is_annotate=="yes"), color="orange", size=3) +
            # add label 
            geom_label_repel(data=subset(dm_annots, is_annotate=="yes"), aes(label=symbol), size=2) +
            labs(x=NULL, y="-log10(FDR pval)") +
            #theme_minimal() +
            theme_bw() +
            theme(
                  legend.position = "none",
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    #axis.title.y = element_markdown(),
                    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
                    panel.border = element_blank()
            ) +
            ggtitle(paste(tr, st, ct, rt, "DM AD genes"))

        manhattandir <- paste0(savedir, "manhattan_plot/")
        filename <- paste(tr, st, ct, rt, "manhattan_plot.png", sep='_')
        savefile <- paste0(manhattandir, filename)
        savefile

        ggsave(savefile, height=6, width=8)
        1
    }

    head(combined_df)
    region_types <- unique(combined_df$region_type)
    cell_types <- unique(combined_df$cell_type)
    summary_types <- unique(combined_df$summary_type)
    traits <- unique(combined_df$trait)

    runs <- expand.grid(region_type=region_types,
                        cell_type=cell_types,
                        summary_type=summary_types,
                        trait=traits
    )

    apply(runs, 1, create_manhattan_plot)

}


plot_results <- function(){
    # create plots to visualize the results
    summary_types <- c('avgs', 'pcs', 'cpgs')
    region_types <- c('full_gene', 'preTSS', 'gene_body', 'promoters')
    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')
    traits <- c('ceradsc', 'cat_braaksc', 'new_apoe', 'new_diag')

    runs <- expand.grid(summary_type=summary_types,
                        region_type=region_types,
                        cell_type=cell_types,
                        trait=traits
    )
    head(runs)
    dim(runs)

    full_array = FALSE

    # load gene symbol file
    datafile <- paste0("/path/to/table/full_biomart_table.csv")
    gene_symbols <- read_csv(datafile) %>%
        separate(Ensembl_ID, c('ensembl_gene_id', 'ensembl_gene_id_version'), sep="\\.")
    head(gene_symbols)

    # load niagds genes
    datafile <- paste0("/path/to/genes/NIAGDS_AD_genes_list.csv")
    niagds_genes <- read_csv(datafile, col_name='Symbol')
    head(niagds_genes)

    # load the dmp results
    row <- runs[1,] # for testing interactively
    load_summary <- function(row, full_array){
        # get parameters for this run
        summary_type <- row[['summary_type']]
        region_type <- row[['region_type']]
        cell_type <- row[['cell_type']]
        trait <- row[['trait']]

        # get the datafile
        datafile <- paste(trait, cell_type, region_type, summary_type, 
                          "unfiltered_limma_results.rds", sep='_'
        )
        if (summary_type == 'cpgs'){
            datafile <- paste(trait, cell_type, region_type, summary_type, 
                              "mapped_unfiltered_limma_results.rds", sep='_'
                              )
        }
        if (summary_type == 'cpgs' & full_array){
            datafile <- paste(trait, cell_type, region_type, summary_type, 
                              "mapped_full_array_unfiltered_limma_results.rds", sep='_'
                              )
        }

        # load the data
        datapath <- paste0(savedir, cleaned_str, datafile)
        datapath
        list.files(savedir)

        if (!file.exists(datapath)){
            print(paste("No datafile for", summary_type, region_type, cell_type, trait))
            return(NA)
        }

        dmp <- readRDS(datapath)
        head(dmp)


        # format the feature column to match across summary types
        if (summary_type == 'cpgs'){
            dmp <- dmp %>%
                mutate(feature = paste0(gene_id, "-", cpg_id))
        } else{
            dmp <- dmp %>%
                rownames_to_column('feature')
        }
        head(dmp)


        # select columns that we care about
        sub_dmp <- dmp %>%
            select(feature, logFC, AveExpr, t, P.Value, bh_pval, bf_pval, B) %>%
            mutate(summary_type = summary_type,
                   region_type = region_type,
                   cell_type = cell_type,
                   full_array=full_array,
                   trait=trait
            )
        head(sub_dmp)

        if (test_bulk){
            sub_dmp <- dmp %>%
                select(feature, logFC, AveExpr, t, P.Value, bh_pval, bf_pval, B, outcome_celltype=cell_type) %>%
                mutate(summary_type = summary_type,
                       region_type = region_type,
                       cell_type = cell_type,
                       full_array=full_array,
                       trait=trait
                )

        }
        head(sub_dmp)

        # make the results nicer
        nice_df <- sub_dmp %>%
            separate(feature, c('gene_id', 'feature_id'), sep='-', fill='right') %>%
            separate(gene_id, c('ensembl_gene_id', 'gene_version'), remove=FALSE, sep="\\.") %>%
            left_join(gene_symbols) %>%
            select(-gene_id) %>%
            rename(gene_id = ensembl_gene_id) %>%
            select(-gene_version, -ensembl_gene_id_version, -full_array) %>%
            select(summary_type:Symbol, feature_id, everything()) %>%
            mutate(niagds_gene = Symbol %in% niagds_genes$Symbol)
        head(nice_df)

        nice_df
    }
    
    # get the results for the matched cpgs
    dmp_res <- apply(runs, 1, load_summary, full_array=FALSE)
    dmp_df <- do.call(rbind, dmp_res)
    head(dmp_df)


    # get the results for the full array cpgs
    unique(dmp_df$summary_type)
    unique(dmp_df$region_type)

    # save the full results
    save_nice_results <- function(region_type){
        runs <- expand.grid(summary_type=summary_types,
                            region_type=region_type,
                            cell_type=cell_types,
                            trait=traits
        )

        print(head(runs))

        # get the results for the matched cpgs
        dmp_res <- apply(runs, 1, load_summary, full_array=FALSE)
        dmp_df <- do.call(rbind, dmp_res)
        head(dmp_df)

        # write the full table of results
        #region_type <- 'full_gene'
        (savefile <- paste0(savedir, region_type, "_full_DM_results.csv"))
        write_csv(dmp_df, savefile)


        # get significant counts
        head(dmp_df)
        sig_df <- dmp_df %>%
            filter(bh_pval < 0.05)
        head(sig_df)

        sig_counts <- sig_df %>%
            rowwise() %>%
            mutate(gene_feature = paste(gene_id, feature_id, sep='-')) %>%
            group_by(summary_type, region_type, cell_type) %>%
            summarize(sig_genes = n_distinct(gene_id),
                      sig_features = n_distinct(gene_feature)
            )
        head(sig_counts)

        # get total counts too
        total_counts <- dmp_df %>%
            rowwise() %>%
            mutate(gene_feature = paste(gene_id, feature_id, sep='-')) %>%
            group_by(summary_type, region_type, cell_type) %>%
            summarize(total_genes = n_distinct(gene_id),
                      total_features = n_distinct(gene_feature)
            )
        head(total_counts)

        # join the sig and total counts together
        all_counts <- total_counts %>%
            left_join(sig_counts) %>%
            mutate(sig_genes = ifelse(is.na(sig_genes), 0, sig_genes)) %>%
            mutate(sig_features = ifelse(is.na(sig_features), 0, sig_features)) 
        head(all_counts)

        (savefile <- paste0(savedir, region_type, "_all_counts.csv"))
        write_csv(all_counts, savefile)
        1
    }
    lapply(region_types, save_nice_results)


    # load files from here to combine together
    (savefiles <- paste0(savedir, region_types, "_all_counts.csv"))
    combined_counts <- lapply(savefiles, read_csv) %>%
        do.call(rbind, .)
    head(combined_counts)
    dim(combined_counts)

    (savefile <- paste0(savedir, "combined_regions_all_counts.csv"))
    write_csv(combined_counts, savefile)

    
    full_df <- NA
    # combine the two results
    combined_df <- rbind(dmp_df, full_df) %>%
        left_join(summary_type_map) %>%
        left_join(cell_type_map) %>%
        left_join(region_type_map)
    head(combined_df)
    dim(combined_df)

    (savefile <- paste0(savedir, "combined_dm_results.rds"))
    saveRDS(combined_df, savefile)

    # types of pvalue adjustment performed
    # plot results for both
    correction_types <- c('bh', 'bf')

    plotting_dir <- paste0(savedir, "plots/")
    dir.create(plotting_dir, showWarnings=FALSE)

    # plot the counts of sig dmps
    correction_type <- correction_types[1]
    plot_dmp_counts_barplot(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=FALSE)
    plot_dmp_counts_barplot(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=TRUE)

    correction_type <- correction_types[2]
    plot_dmp_counts_barplot(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=FALSE)
    plot_dmp_counts_barplot(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=TRUE)

    
    # plot the counts of sig genes
    correction_type <- correction_types[1]
    plot_gene_counts_barplot(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=FALSE)
    plot_gene_counts_barplot(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=TRUE)

    correction_type <- correction_types[2]
    plot_gene_counts_barplot(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=FALSE)
    plot_gene_counts_barplot(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=TRUE)


    # check pcs
    correction_type <- correction_types[1]
    check_pcs(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=TRUE)

    # plot manhattan plot
    correction_type <- correction_types[1]
    plot_manhattan(combined_df, correction_type, summary_type_colors, plotting_dir, full_array=FALSE)

    head(combined_df)
}


main <- function(){
    # define parameter choices for all runs
    region_types <- c('gene_body', 'preTSS', 'full_gene', 'promoters')
    summary_types <- c('pcs', 'avgs', 'cpgs')
    traits <- c('ceradsc', 'cat_braaksc', 'new_apoe', 'new_diag')
    cleaned_strs <- c("cleaned_global_")
    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')

    runs <- expand.grid(region_type=region_types,
                        summary_type=summary_types,
                        trait=traits,
                        cleaned_str=cleaned_strs,
                        cell_type=cell_types
    ) %>%
        mutate(rn = row_number())
    head(runs, 10)
    
    row <- runs[1,] # for testing interactively
    process_row <- function(row){
        # get the parameters for this run
        region_type <- row[['region_type']] %>% as.character()
        summary_type <- row[['summary_type']] %>% as.character()
        cleaned_str <- row[['cleaned_str']] %>% as.character()
        trait <- row[['trait']] %>% as.character()
        cell_type <- row[['cell_type']] %>% as.character()
        rownum <- row[['rn']] %>% as.character()

        print(paste("Processing", region_type, summary_type, cleaned_str, trait, cell_type))
        print(paste("Row", rownum, "out of", nrow(runs)))

        # load phenotype data
        pheno <- load_pheno(cell_type)
        head(pheno)


        # load methylation data
        meth <- load_meth(region_type, summary_type, cleaned_str, cell_type)
        head(meth)
        head(meth)[1:10]


        # reduce the cpgs to only those used to make the averages and regional pcs
        if (summary_type == 'cpgs'){
            print("Matching cpgs to region annotations")

            # match the cpgs to the regions
            # load cpg-gene map for region
            datafile <- paste0(datadir,
                               "/", region_type, "_", cell_type, "_gene_cpg_map.rds")
            cpg_map <- readRDS(datafile)
            head(cpg_map)

            # get cpg positions and id
            datafile <- paste0(datadir,
                               cell_type, "_cpg_position_map.rds")
            pos_map <- readRDS(datafile) %>%
                rename(cpg_id = cpg)
            head(pos_map)

            # combine annots
            ext_map <- cpg_map %>%
                left_join(pos_map)
            head(ext_map)

            head(meth)

            keep_probes <- intersect(rownames(meth), ext_map$cpg_id)
            length(keep_probes)

            sub_meth <- meth[keep_probes,]
            dim(meth)
            dim(sub_meth)

            meth <- sub_meth
        }


        # get global PCs to use as covariates
        global_pcs = NA
        if (include_global){
            # select the appropriate global pcs
            if (cell_type != 'bulk'){
                pcs_file <- paste0(svddir, cleaned_str, "cleaned_", cell_type, "_global_pcs.rds")
            } else{
                pcs_file <- paste0(svddir, "cleaned_bulk_global_pcs.rds")
            }
            pcs_file

            # read in the global pcs that were computed in the 
            # 02_batch_correction script in the 02_deconvolution directory
            pcs <- readRDS(pcs_file)
            pcs <- pcs[-1,]
            head(pcs)

            if (protect_global){
                # do not include PCs if they have no significant correlation with
                # the outcome trait, batch, pmi, and study
                print(paste("Protecting global PCs"))

                 #filter pcs to only those that we want to remove
                 #load lm pvals for global PCs
                savefile <- paste0(svddir, cleaned_str, "cleaned_", cell_type, "_pc_lm_pvals.rds")
                lm_pval <- readRDS(savefile)
                head(lm_pval)

                # find the pcs to protect
                thresh <- 0.1
                protect_pcs <- lm_pval %>%
                    select(-pct_var, -pval) %>%
                    spread(trait, bf_pval)

                protect_pcs['trait'] <- protect_pcs[trait]
                protect_pcs <- protect_pcs %>%
                    select(pc, trait, meth.batch, pmi, Study) %>%
                    filter(
                           trait < thresh,
                           meth.batch >= thresh,
                           pmi >= thresh,
                           Study >= thresh
                    ) %>%
                    mutate(pc = str_replace(pc, 'pc_', 'PC'))
                head(protect_pcs)

                # these pcs should be removed from covariates
                head(pcs)
                remove_pcs <- setdiff(colnames(pcs), protect_pcs$pc)
                remove_pcs

                global_pcs <- pcs[,remove_pcs] 
            } else{
                print(paste("NOT protecting global PCs"))
                global_pcs <- pcs
            } 
        }
        head(global_pcs)

        # run dmp analysis
        dmp_res <- run_dmp(meth, pheno, global_pcs, region_type, summary_type, cleaned_str, trait, cell_type, shuffle_pheno, shuffle_num=1)
        print("-----------------------------------------------")

        1
    }

    dmp_res <- apply(runs, 1, process_row)

    # check the results
    check_results()

    # plot the final results
    plot_results()
}

main()
