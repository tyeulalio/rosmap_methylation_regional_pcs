library(wateRmelon)
library(EpiSCORE)
library(TCA)
library(RNOmni)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)

# Use EpiSCORE to estimate cell type proportions for ROSMAP methylation data
# Use Tensor Composition Analysis to deconvolve ROSMAP methylation data
# Run this on betas from the normalized and QC'ed beta values

savedir <- "/path/to/savedir"

load_methylation <- function(thresh){
    # thresh = allele frequency threshold used in preprocessing

    # preprocessed ROSMAP methylation data
    datafile <- paste0("../../output/01_preprocessing/methylation/01_preprocessing_methylation/snp", thresh, "_window0_filtered_betas.rds")

    betas <- readRDS(datafile)
    head(betas)
    dim(betas)

    betas
}

load_clinical <- function(){
    # load phenotype data
    datafile <- paste0("../../output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds")
    targets <- readRDS(datafile)
    head(targets)

    # grab mapping from meth colname and specimenID
    datafile <- paste0("../../data/ROSMAP_data/Epigenetics/Epigenetics (DNA methylation array)/SYNAPSE_METADATA_MANIFEST.tsv")
    manifest <- read_tsv(datafile)
    head(manifest)

    # map colname to specimenID
    name_map <- manifest %>%
        filter(fileFormat == 'idat') %>%
        select(specimenID, name) %>%
        mutate(name = gsub("_[A-Za-z]{3}\\.idat", "", name)) %>%
        unique() %>%
        filter(!is.na(specimenID))
    head(name_map)

    # map specimen to individualID
    specimen_map <- targets %>%
        select(individualID, specimenID=meth.specimenID) %>%
        inner_join(name_map) %>%
        rename(meth.specimenID=specimenID)
    head(specimen_map)

    # connect to clincal
    mapped_clinical <- targets %>%
        left_join(specimen_map)
    head(mapped_clinical)

    savefile <- paste0(savedir, "mapped_clinical.rds")
    savefile
    saveRDS(mapped_clinical, savefile)

    dim(mapped_clinical)

    mapped_clinical
}

format_data <- function(betas, mapped_clinical){
    # format the betas data and clinical data
    head(betas)
    dim(betas)
    head(mapped_clinical)

    # put both data in same order
    head(betas)
    head(mapped_clinical)

    # get the names from both dataframes
    overlapped_names <- intersect(colnames(betas), mapped_clinical$name) %>% na.omit()
    length(overlapped_names)

    # set both dataframes in the same order
    ordered_clinical <- mapped_clinical[match(overlapped_names, mapped_clinical$name),]
    ordered_betas <- betas[,overlapped_names]


    if (!identical(colnames(ordered_betas), ordered_clinical$name)) stop("ERROR: ORDERING DOES NOT MATCH")
    

    # set the formatted dataframe and set rownames as indivudalID
    formatted_betas <- ordered_betas
    colnames(formatted_betas) <- ordered_clinical$individualID
    head(formatted_betas)

    head(ordered_clinical)
    names(ordered_clinical)


    # determine relationship between covariates
    # separate categorical and numerical variables
    factor_covariates <- c('sex', 'race', 'spanish', 'cogdx', 'diagnosis', 'apoe4', 'meth.batch',
                           'meth.platform', 'meth.Sentrix_ID', 'meth.Sentrix_Row_Column', 
                           'meth.Sample_Plate', 'meth.Sample_Well')
    cont_covariates <- c('age_death', 'pmi', 'educ')

    # remove +'s from variables - mainly for age variable
    ordered_clinical <- data.frame(lapply(ordered_clinical,
                                          function(x){x <- sapply(x,function(y){str_replace_all(as.character(y),'\\+','')})}))
    rownames(ordered_clinical) <- ordered_clinical$individualID
    head(ordered_clinical)

    # Convert factor covariates to factors and numeric covariates to numeric
    ordered_clinical[,factor_covariates] = lapply(ordered_clinical[,factor_covariates], factor)
    ordered_clinical[,cont_covariates] = lapply(ordered_clinical[,cont_covariates], as.character)
    ordered_clinical[,cont_covariates] = lapply(ordered_clinical[,cont_covariates], as.numeric)

    head(ordered_clinical)

    formatted_clinical <- ordered_clinical

    # save these formatted betas
    savefile <- paste0(savedir, "formatted_betas.rds")
    savefile
    saveRDS(formatted_betas, savefile)

    formatted_dat <- list(formatted_betas=formatted_betas, formatted_clinical=formatted_clinical)
    formatted_dat
}


get_cell_proportions <- function(formatted_betas, formatted_clinical, thresh){
    # get the proportions of cell types in the data
    head(formatted_betas)

    # use the noob preprocessed data
    load_noob <- function(){
        datafile <- paste0("../..//output/01_preprocessing/methylation/01_preprocessing_methylation/noob_preprocessed_meth.rds")
        datafile
        noob_meth <- readRDS(datafile)
        head(noob_meth)

        # adjust the sample names
        head(formatted_clinical)
        ordered_samples <- intersect(colnames(noob_meth), formatted_clinical$name)
        length(ordered_samples)

        # match the order of betas and clinical
        ordered_meth <- noob_meth[,ordered_samples]
        ordered_clinical <- formatted_clinical[match(ordered_samples, formatted_clinical$name),]

        # make sure order is the same
        stopifnot(identical(colnames(ordered_meth), ordered_clinical$name))

        colnames(ordered_meth) <- ordered_clinical$individualID
        head(ordered_meth)

        ordered_meth
    }
    ordered_meth <- load_noob()
    head(ordered_meth)


    # load reference matrix for brain from EpiSCORE
    data(BrainRef) # should load exprefBrain.m and mrefBrain.m
    head(mrefBrain.m)
    head(exprefBrain.m)

    brain_ref <- mrefBrain.m
    head(brain_ref)

    # summarise betas 
    # computes average beta across cpgs 200bp upstream of TSS
    # cols = samples
    # rows = cpgs

    head(ordered_meth )
    # try remove low variance probes here
    sds <- rowSds(as.matrix(ordered_meth))
    summary(sds)

    # set the threshold here
    sd_thresh = 0
    table(sds >= sd_thresh)

    var_betas <- ordered_meth[sds >= sd_thresh,]

    betas_mat <- as.matrix(var_betas)
    dim(betas_mat)
    head(betas_mat)

    # adjust array type here
    arraytype = "850k"

    # get averages for each reference gene
    avg_betas = constAvBetaTSS(betas_mat, type=arraytype)
    head(avg_betas)


    # decompose data using EpiSCORE (estimate cell type proportions)
    head(brain_ref)
    celltype_proportions <- wRPC(avg_betas, brain_ref, maxit=500)
    
    names(celltype_proportions)
    head(celltype_proportions$estF)
    head(celltype_proportions$ref)

    # save the data
    savefile <- paste0(savedir, precleaned_str, "snpThresh", thresh, "_episcore_estimated_celltypes.rds")
    savefile 
    saveRDS(celltype_proportions$estF, savefile)

    savefile <- paste0(savedir, precleaned_str, "snpThresh", thresh, "_episcore_reference_panel.rds")
    savefile 
    saveRDS(celltype_proportions$ref, savefile)


    celltype_proportions$estF
}

deconvolve_data <- function(formatted_betas, proportions_df, formatted_clinical){
    # deconvolve data using Tensor Composition Analysis
    head(proportions_df)
    summary(proportions_df$Oligo) # very small, combine with oligo
    summary(proportions_df$OPC) # very small, combine with oligo
    summary(proportions_df$Microglia) # zero microliga -- remove

    # sum the proportions of oligodendorcytes and OPCs
    proportions_df <- proportions_df %>%
        mutate(oligo_opc = Oligo + OPC)
    head(proportions_df)

    proportions_df %>%
        filter(oligo_opc != OPC) %>%
        head() # only 4 samples with Oligo

    # select the proportions that we'll use for deconvolution
    sub_props <- proportions_df %>%
        select(neuron, oligo_opc, astro=Astro, endo=Endo)
    head(sub_props)

    savefile <- paste0(savedir, precleaned_str, "formatted_celltype_proportions.rds")
    saveRDS(sub_props, savefile)

    head(formatted_betas)
    head(sub_props)



    # decide whether to add batch in the tca model
    # c1 = covariates taht affect methylation at the cell type level (age, gender)
    # c2 = covariates that are not expected to affect meth at the celll type level,
    # but rather may affect tissue-level meth (batch, technical variation)
    head(formatted_clinical)
    identical(colnames(formatted_betas), formatted_clinical$individualID)

    sub_clinical <- formatted_clinical %>%
        filter(!is.na(age_death)) %>%
        mutate(study_num = as.numeric(as.factor(Study)))

    table(sub_clinical$study_num)

    c2_vars = c('meth.batch', 'pmi', 'study_num')
        
    # set up c2 design matrix
    # replace missing pmi with the mean
    sub_clinical <- sub_clinical %>%
        mutate(pmi = ifelse(is.na(pmi), mean(pmi, na.rm=TRUE), pmi))

    f = paste0("~ ", paste(c2_vars, collapse = ' + ')) %>% as.formula()
    f
    c2_design <- model.matrix(f, data=sub_clinical)[,-1]
    head(c2_design)

    rownames(c2_design) = sub_clinical$individualID
    colnames(c2_design)

    head(c2_design)
    dim(c2_design)


    # select betas and props base don indivudals
    sub_betas = formatted_betas[,sub_clinical$individualID]
    dim(sub_betas)
    
    sub_props <- sub_props[sub_clinical$individualID,]
    head(sub_props)



    # run the tca deconvolution
    # betas - cols=samples, rows=cpgs
    # proportions = cols=celltypes, rows=samples
    tca_out <- tca(X=sub_betas,
                   W=sub_props,
                   C2=c2_design,
                   #parallel = TRUE, num_cores = 16
                   constrain_mu=TRUE,
                   refit_W=TRUE
    )

    # save the results
    savefile <- paste0(savedir, precleaned_str, "tca_out.rds")
    savefile
    saveRDS(tca_out, savefile)


    # project from original mixed betas matrix to celltype-specific tensor
    # this step uses multilinear algebra
    head(sub_betas)
    betas_mat <- as.matrix(sub_betas)
    tca_tensor <- tensor(betas_mat, tca_out, parallel=FALSE)

    savefile <- paste0(savedir, precleaned_str, "tca_tensor.rds")
    savefile
    saveRDS(tca_tensor, savefile)

    1
}

main <- function(){
    # run through the main workflow
    print("Loading data")

    # load normalized methylation
    thresh <- 0
    betas <- load_methylation(thresh)
    dim(betas)

    # load the phenotype data
    mapped_clinical <- load_clinical()

    # format the data 
    formatted_dat <- format_data(betas, mapped_clinical)
    names(formatted_dat)
    formatted_betas <- formatted_dat$formatted_betas
    formatted_clinical <- formatted_dat$formatted_clinical

    head(formatted_clinical)
    dim(formatted_betas)

    savefile <- paste0(savedir, "formatted_clinical.rds")
    saveRDS(formatted_clinical, savefile)

    formatted_clinical <- readRDS(savefile)


    # get the cell type proportions 
    print(paste("Getting cell type proportions", Sys.time()))
    celltype_proportions <- get_cell_proportions(formatted_betas, formatted_clinical, thresh)


    # deconvolve the data
    print(paste("Deconvolving data", Sys.time()))
    proportions_df <- as.data.frame(celltype_proportions)
    head(proportions_df)

    deconvolve_data(formatted_betas, proportions_df, formatted_clinical)

}

main()
