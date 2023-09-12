library(isva) # required for smartSVq
library(wateRmelon)
library(RNOmni)
library(PCAtools)
library(limma)
library(ggplot2)
library(tidyverse)


# This script removes effects from confouding variables from the methylation data
# Run this on mvals from the normalized and QC'ed beta values

datadir <- "/path/to/data"
savedir <- "/path/to/save/directory"

load_methylation <- function(cell_type){
    # load the methylation data

    # for bulk data, just use the preprocessed methylation
    if (cell_type == 'bulk'){
        datafile <- paste0(datadir, "formatted_betas.rds")

        betas <- readRDS(datafile)
        head(betas)
        return(betas)
    }

    # for other cell types, get the deconvolved beta values
    datafile <- paste0(datadir, "tca_tensor.rds")
    datafile
    tca_tensor <- readRDS(datafile)

    # get the ordered celltypes 
    datafile <- paste0(datadir, "tca_out.rds")
    tca_out <- readRDS(datafile)

    # set the cell type names
    names(tca_tensor) <- colnames(tca_out$W)

    # grab the cell type betas
    betas <- tca_tensor[[cell_type]]
    head(betas)[1:10]

    betas
}

correlate_traits <- function(clinical_cts){ 
    # use this to correlate phenotypes 
    # with one another to see how they're related
    # the main script doesn't call this but run this
    # interactively to get the results if needed
    head(clinical_cts)
    names(clinical_cts)

    # select phenotypes that we care about
    # format them as numeric or factors
    sub_clinical <- clinical_cts %>%
        select(individualID, study=Study, msex, educ, apoe_genotype, apoe4, age_death, pmi, braaksc, ceradsc, cogdx, dcfdx_lv,
               astro, endo, neuron, oligo_opc, new_diag, new_apoe) %>%
        # remove + from age_death and make numeric
        mutate(age_death = as.numeric(str_remove(age_death, "\\+"))) %>%
        # create categorical traits
        mutate(braaksc_cat = as.numeric(braaksc >= 5)) %>%
        mutate(ceradsc_cat = as.factor(as.character(ceradsc))) %>%
        mutate(dcfdx_lv_cat = as.factor(as.character(dcfdx_lv))) %>%
        mutate(cogdx_cat = as.factor(as.character(cogdx))) %>%
        # apoe should not be continuous
        mutate(apoe_genotype = as.factor(as.character(apoe_genotype)),
               apoe4 = as.factor(as.character(apoe4))
        ) %>%
        column_to_rownames('individualID')
    head(sub_clinical)
    str(sub_clinical)

    # need to run correlations between each of the traits
    (traits <- names(sub_clinical))

    # create pairwise comparisons
    runs <- expand.grid(trait1=traits, trait2=traits)
    head(runs)
    runs
    
    # for testing
    row <- runs[6,]
    trait1 <- "ceradsc_cat"
    trait2 <- "ceradsc_cat"

    get_correlation <- function(row){
        # correlate pairs of traits
        trait1 <- row[['trait1']] %>% as.character()
        trait2 <- row[['trait2']] %>% as.character()

        #print(paste("Correlating", trait1, trait2))

        # get the classes of each trait
        (class1 <- class(sub_clinical[[trait1]]))
        (class2 <- class(sub_clinical[[trait2]]))

        #print(paste(class1, class2))

        isfactor1 <- class1 == "factor" | class1 == "character"
        isfactor2 <- class2 == "factor" | class2 == "character"

        # case 1: both numeric
        if (!(isfactor1 | isfactor2)){ 
            res <- cor.test(sub_clinical[[trait1]], sub_clinical[[trait2]], method='spearman')
            pval <- res$p.value
            method <- "spearman"
        }

        # case 2: both character
        if (isfactor1 & isfactor2){ 
            res <- chisq.test(sub_clinical[[trait1]], sub_clinical[[trait2]])             
            pval <- res$p.value
            method <- "chisq"
        }

        # case 3/4: one numeric, one character
        if (isfactor1 & !isfactor2){
           res <- kruskal.test(sub_clinical[[trait2]], sub_clinical[[trait1]]) 
            pval <- res$p.value
            method <- "kruskal-wallis"
        }
        if (!isfactor1 & isfactor2){
           res <- kruskal.test(sub_clinical[[trait1]], sub_clinical[[trait2]]) 
            pval <- res$p.value
            method <- "kruskal-wallis"
        }
        
        tibble(trait1=trait1,
               trait2=trait2,
               pval=pval,
               method=method
        )
    }
    corrs <- apply(runs, 1, get_correlation)
    corrs_df <- do.call(rbind, corrs)
    head(corrs_df)

    # add adjusted p-values
    formatted_corrs <- corrs_df %>%
        mutate(bf_pval = p.adjust(pval, method='bonferroni'))
    head(formatted_corrs)

    # save results
    savefile <- paste0(savedir, "trait_trait_correlations.rds")
    saveRDS(formatted_corrs, savefile)

    # plot the results as heatmap
    plotdir <- paste0(savedir, "plots/")
    (savefile <- paste0(plotdir, "trait_correlation_heatmap.png"))
    traits
    ordered_traits <- c('age_death', 'msex', 'educ', 'study', 'pmi', 'astro', 'endo',
                        'neuron', 'oligo_opc', 'apoe_genotype', 'apoe4', 'new_apoe',
                        'cogdx', 'cogdx_cat', 'dcfdx_lv', 'dcfdx_lv_cat', 'new_diag',
                        'braaksc', 'braaksc_cat', 'ceradsc', 'ceradsc_cat'
    )
    p <- formatted_corrs %>%
        mutate(trait1 = factor(trait1, levels=ordered_traits)) %>%
        mutate(trait2 = factor(trait2, levels=ordered_traits)) %>%
        ggplot(aes(x=trait1, y=trait2)) +
        geom_tile(aes(fill=bf_pval)) +
        geom_text(aes(label=round(bf_pval,2))) +
        scale_fill_continuous(low='red', high='white', limits=c(0,0.1)) +
        theme_bw() +
        theme(text=element_text(size=20),
              axis.text.x=element_text(angle=30, vjust=1, hjust=1) )
    savefile
    ggsave(p, file=savefile, width=15, height=8)

}


load_clinical <- function(){
    # grab clinical data
    datafile <- paste0("../../output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds")
    clinical <- readRDS(datafile)
    head(clinical)

    # add cell type proportion data
    datafile <- paste0(datadir, "tca_out.rds")
    tca_out <- readRDS(datafile)
    names(tca_out)

    props <- tca_out$W
    head(props)

    # create data frame
    celltype_props <- props %>%
        as.data.frame() %>%
        rownames_to_column('individualID')
    head(celltype_props)

    # join proportions to phenotype traits
    clinical_cts <- clinical %>%
        left_join(celltype_props) 
    names(clinical_cts)
    head(clinical_cts)

    # create new clinical groups
    new_clinical <- clinical_cts %>%
        mutate(new_diag = ifelse(cogdx %in% 2:5 & braaksc >= 4 & ceradsc <= 2, "AD", "Other"),
               new_diag = ifelse(cogdx == 1 & braaksc <=3 & ceradsc >= 3, "Control", new_diag)) %>%
        mutate(new_apoe = ifelse(apoe_genotype %in% c(33,34), apoe_genotype, NA)) %>%
        mutate(new_braaksc = ifelse(braaksc == 6, NA, braaksc))
    head(new_clinical)
    
    clinical_cts <- new_clinical

    return(clinical_cts)

    # create plot for correlation of traits
    # only need to do this once 
    correlate_traits(clinical_cts)
}

format_data <- function(betas, mapped_clinical, cell_type){
    # format the betas data and clinical data
    head(betas)
    head(mapped_clinical)

    # put both data in same order
    overlapped_names <- intersect(colnames(betas), mapped_clinical$individualID) %>% na.omit()
    length(overlapped_names)

    ordered_clinical <- mapped_clinical[match(overlapped_names, mapped_clinical$individualID),]
    ordered_betas <- betas[,overlapped_names]

    if (!identical(colnames(ordered_betas), ordered_clinical$individualID)) stop("ERROR: ORDERING DOES NOT MATCH")
    

    formatted_betas <- ordered_betas
    colnames(formatted_betas) <- ordered_clinical$individualID
    head(formatted_betas)

    head(ordered_clinical)
    names(ordered_clinical)


    # determine relationship between covariates
    factor_covariates <- c('sex', 'race', 'spanish', 'cogdx', 'diagnosis', 'apoe4', 'meth.batch', 'apoe_genotype',
                           'meth.platform', 'meth.Sentrix_ID', 'meth.Sentrix_Row_Column', 'meth.Sample_Plate', 'meth.Sample_Well')
    cont_covariates <- c('age_death', 'pmi', 'educ', 'neuron', 'astro', 'endo', 'oligo_opc')

    ordered_clinical <- data.frame(lapply(ordered_clinical,function(x){x <- sapply(x,function(y){str_replace_all(as.character(y),'\\+','')})}))
    rownames(ordered_clinical) <- ordered_clinical$individualID
    head(ordered_clinical)

    # Convert factor covariates to factors
    ordered_clinical[,factor_covariates] = lapply(ordered_clinical[,factor_covariates], factor)
    ordered_clinical[,cont_covariates] = lapply(ordered_clinical[,cont_covariates], as.character)
    ordered_clinical[,cont_covariates] = lapply(ordered_clinical[,cont_covariates], as.numeric)

    head(ordered_clinical)

    # create categorical braak based on braak score
    ordered_clinical <- ordered_clinical %>%
        mutate(cat_braaksc=braaksc >=5) %>%
        mutate(cat_braaksc=as.numeric(cat_braaksc)) %>%
        mutate(bin_apoe4=as.numeric(apoe_genotype %in% c('24', '34', '44')))
    table(ordered_clinical$cat_braaksc, ordered_clinical$braaksc)

    formatted_clinical <- ordered_clinical

    formatted_dat <- list(formatted_betas=formatted_betas, formatted_clinical=formatted_clinical)
    formatted_dat
}

process_betas <- function(formatted_betas, cell_type){
    # process the beta values
    # - remove low variance cpgs
    # - convert to mvals
    # - center and scale cpgs

    # remove low variance cpgs
    # this may have already been applied in earlier script but just do it again
    # this is similar to what execute.variability.removal from the rnbeads pacakge
    var_betas <- formatted_betas[apply(formatted_betas, 1, var, na.rm=TRUE) != 0,]
    dim(var_betas)
    dim(formatted_betas)
    head(var_betas)

    nrow(formatted_betas) - nrow(var_betas) # number of cpgs lost to low variance
    sum(is.na(var_betas))

    # return the appropriate meth form
    if (meth_conversion == 'INT'){
        print("Using INT conversion")

        # converts methylation using inverse normal transform
        # inverse normal transform
        int_meth <- apply(var_betas, 1, RankNorm) %>% 
            t() %>%
            as.data.frame()
        head(int_meth)

        return(int_meth)
    }
    if (meth_conversion == 'mvals'){
        print("Using mval conversion")
        # transform methylation to m-values
        print("scaling betas after TCA")
        betas_matrix <- as.matrix(var_betas)
        betas_matrix[betas_matrix < 0] <- 0
        betas_matrix[betas_matrix > 1] <- 1
        betas_matrix <- (betas_matrix * (ncol(betas_matrix) - 1) + 0.5) / ncol(betas_matrix)
        sum(is.na(betas_matrix))

        # convert to mvals
        mvals <- beta2m(betas_matrix) %>%
            as.data.frame()
        head(mvals)
        sum(is.na(mvals))

        return(mvals)
    }
    
    1
}


run_pca <- function(mvals, cell_type, cv_savedir, datatype=""){
    # perform PCA to get PCs without protecting any traits
    head(mvals)

    # get PCs
    # rows = variables
    # cols = samples
    pca_res <- PCAtools::pca(mvals)

    # get dim using gavish donoho
    # rows = variables
    # cols = observations
    gv_res <- chooseGavishDonoho(mvals, var.explained=pca_res$sdev^2, noise=1)
    gv_res

    mp_res <- chooseMarchenkoPastur(mvals, var.explained=pca_res$sdev^2, noise=1)
    mp_res

    print(paste("gavish donoho dim:", gv_res[[1]], "marchenko pastur dim:", mp_res[[1]]))

    # get the pcs
    names(pca_res)
    pcs <- pca_res$rotated
    head(pcs)

    # get the estimated dimension
    est_dim <- gv_res[[1]]
    est_dim

    # subset PCs to the estimated dimension 
    sig_pcs <- pcs[,1:est_dim,drop=FALSE] %>%
        as.data.frame()
    head(sig_pcs)
    dim(sig_pcs)

    # get pct var explained
    eig_sq <- pca_res$sdev ^ 2
    pct_var <- eig_sq / sum(eig_sq)
    subset_pct_var <- pct_var[1:est_dim]

    # store the percent variance explained together with the PCs 
    formatted_pcs <- rbind(subset_pct_var, sig_pcs)
    rownames(formatted_pcs)[1] <- "percent_variance"
    head(formatted_pcs)

    # print just to see what it is
    mp_var <- sum(pct_var[1:mp_res[[1]]])
    print(paste("Variance explained, GV:", round(sum(subset_pct_var), 2), "MP:", round(sum(mp_var), 2)))

    # save the results
    print("Saving files")
    savefile <- paste0(cv_savedir, datatype, cell_type, "_global_pcs.rds")
    savefile
    saveRDS(formatted_pcs, savefile)

    1
}


correlate_pcs <- function(formatted_clinical, cell_type, cv_savedir, datatype=""){
    # correlate PCs with trait
    targets <- formatted_clinical
    head(targets)
    dim(targets)


    # traits that we'll test against 
    names(targets)
    traits <- c('individualID', 'Study', 'educ', 'age_at_visit_max', 'age_first_ad_dx', 
                'age_death', 'cts_mmse30_first_ad_dx', 'cts_mmse30_lv', 'pmi', 
                'meth.batch', 'sex', 'apoe4', 'apoe_genotype', 'bin_apoe4',
                'cat_braaksc', 'braaksc', 'ceradsc',
                'cogdx', 'dcfdx_lv', 'diagnosis', 'new_diag', 'new_apoe',
                'neuron', 'oligo_opc', 'astro', 'endo'
    )

    # subset the targets
    target_traits <- targets[traits]
    head(target_traits)


    # read in the needed files
    pc_file <- paste0(cv_savedir, precleaned_str, datatype, cell_type, "_global_pcs.rds")
    pcs <- readRDS(pc_file)
    head(pcs)

    # get the percetn variance explained
    pct_var <- pcs[1,]

    # get overlap of inividuals from pcs and traits dataframes
    ordered_samples <- intersect(rownames(pcs), target_traits$individualID)

    # make sure samples are ordered appropriately    
    ordered_traits <- target_traits[match(ordered_samples, target_traits$individualID),]
    head(ordered_traits)


    # put pcs in dataframe
    pc_df <- pcs[ordered_samples,] %>%
        as.data.frame()

    # store the number of pcs
    n.pc <- ncol(pcs)

    # store number of traits
    n.traits <- length(traits) - 1

    # matrix to store correlation p-values
    pvals_mat <- matrix(nrow=n.pc, ncol=n.traits)
    cols <- rep(0,n.traits)

    # correlte the PCs with each of the traits
    pc_num <- 1
    trait_num <- 1
    for (pc_num in 1:n.pc){
        for (trait_num in 1:n.traits){
            trait <- traits[trait_num+1]

            # kruskal wallis test for categorical traits
            # Spearman test for numeric traits
            if (class(ordered_traits[,trait]) != 'numeric'){
                pvals_mat[pc_num, trait_num] <- kruskal.test(pc_df[,pc_num] ~ ordered_traits[,trait])$p.value
            } else{
                pvals_mat[pc_num, trait_num] <- cor.test(pc_df[,pc_num], ordered_traits[,trait], method='spearman',
                                                         exact=FALSE)$p.value
            }

            cols[trait_num] <- trait
        }
    }
    colnames(pvals_mat) <- cols
    pvals_mat

    # check correltion with controls
    head(pc_df)
    head(ordered_traits)


    # store a corrected and uncorrected version of pvals
    pvals_long <- pvals_mat %>%
        as.data.frame() %>%
        mutate(pc=paste("pc", row_number(), sep='_')) %>%
        gather(key='trait', value='pval', -pc) %>%
        group_by(trait) %>%
        mutate(bf_pval = p.adjust(pval, method='bonferroni')) %>%
        mutate(pct_var=t(pct_var)[,1])
    head(pvals_long)


    # save the output
    savefile <- paste0(cv_savedir, precleaned_str, datatype, cell_type, "_pc_lm_pvals.rds")
    if (clean_global & datatype == 'cleaned_'){
        savefile <- paste0(cv_savedir, cleaned_str, precleaned_str, datatype, cell_type, "_pc_lm_pvals.rds")
    }
    savefile
    saveRDS(pvals_long, savefile)

    1
}

limma_remove_batch <- function(formatted_mvals, formatted_clinical, cell_type, cv_savedir){
    # move traits that are highly correlated with outcome
    head(formatted_clinical)
    head(formatted_mvals)

    cleaned_mvals <- NA
    # variables to remove
    remove_vars <- c('batch', 'pmi', 'study_num', 'msex')
    print(paste(c("removing", remove_vars), collapse=' '))

    # remove global PCs too
    global_pcs <- paste0("PC", 1:4)
    remove_vars <- c(remove_vars, global_pcs)
    remove_vars

    # load global PCs
    datafile <- paste0(savedir, cell_type, "_global_pcs.rds")
    global_pcs <- readRDS(datafile)[-1,] %>%
        rownames_to_column('individualID')
    head(global_pcs)

    # covert study to numeric
    temp <- formatted_clinical %>%
        mutate(study_num = as.numeric(as.factor(Study))) %>%
        left_join(global_pcs)
    head(temp)
    formatted_clinical <- temp

    # get ordered samples
    ordered_samples <- colnames(formatted_mvals)
    ordered_samples

    # make sure the batch variable is ordered correctly
    ordered_batch <- formatted_clinical 
    head(ordered_batch)
    rownames(ordered_batch) <- NULL
    ordered_batch <- ordered_batch %>%
        rename(batch=meth.batch) %>%
        select(individualID, all_of(remove_vars)) %>% 
        mutate(pmi = ifelse(is.na(pmi), mean(pmi, na.rm=TRUE), pmi)) %>%
        column_to_rownames('individualID')
    ordered_batch <- ordered_batch[ordered_samples,remove_vars]
    head(ordered_batch)

    stopifnot(identical(rownames(ordered_batch), colnames(formatted_mvals)))

    # separate batch from the other variables
    covariates_matrix = ordered_batch %>%
        select(-batch) %>%
        data.matrix()
    head(covariates_matrix)

     #remove batch variable
     #rows = probes, cols = samples
    cleaned_mvals <- removeBatchEffect(x=formatted_mvals, 
                                       batch=ordered_batch$batch, 
                                       covariates=covariates_matrix
    ) 
    head(cleaned_mvals)

    # save this
    savefile <- paste0(cv_savedir, cleaned_str, cell_type, "_cleaned_meth.rds")
    savefile
    saveRDS(cleaned_mvals, savefile)

    cleaned_mvals
}

main <- function(){
    # process each cell type separately
    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')

    cell_type <- cell_types[celltype_idx]
    process_celltype <- function(cell_type){
        print(paste("Processing", cell_type, Sys.time()))

        print("Loading data")
        formatted_mvals = NA
        formatted_clinical = NA
        savefile <- paste0(savedir, precleaned_str, cell_type, "_formatted_meth.rds")
        # load the cpg methylation data
        betas <- load_methylation(cell_type)
        sum(is.na(betas))

        print("TCA output betas dim")
        print(paste(dim(betas)))

        # load the phenotype data
        mapped_clinical <- load_clinical()
        head(mapped_clinical)

        # format the data
        formatted_dat <- format_data(betas, mapped_clinical, cell_type)
        names(formatted_dat)
        formatted_betas <- formatted_dat$formatted_betas
        formatted_clinical <- formatted_dat$formatted_clinical
        sum(is.na(formatted_betas))

        # convert betas to mvals

        # select the type of methylation conversion
        formatted_mvals <- process_betas(formatted_betas, cell_type)
        sum(is.na(formatted_mvals))

        print("Processed betas")
        print(paste(dim(formatted_mvals)))

        savefile <- paste0(savedir, precleaned_str, cell_type, "_formatted_meth.rds")
        savefile
        saveRDS(formatted_mvals, savefile)

        savefile <- paste0(savedir, cell_type, "_formatted_clinical.rds")
        savefile
        saveRDS(formatted_clinical, savefile)

        head(formatted_mvals)
        dim(formatted_mvals)

        cv_savedir <- savedir
        cv_savedir
        dir.create(cv_savedir)

        # run svd on the mvals
        print("Running PCA")
        datatype = ""
        run_pca(formatted_mvals, cell_type, cv_savedir, datatype)

        print("Correlate PCs with trait")
        correlate_pcs(formatted_clinical, cell_type, cv_savedir, datatype)


        if (cleaning){
            print("Removing batch")
            cleaned_mvals <- limma_remove_batch(formatted_mvals, formatted_clinical, cell_type, cv_savedir)

             #run svd again for plotting
            datatype="cleaned_"
            print("Running SVD on cleaned data")
            run_pca(cleaned_mvals, cell_type, cv_savedir, datatype="cleaned_")

            print("Correlate PCs with trait")
            correlate_pcs(formatted_clinical, cell_type, cv_savedir, datatype="cleaned_")
        }

        1
    }

    cleaning=TRUE
    lapply(cell_types, process_celltype)

    process_celltype(cell_type)

}

main()

