library(RnBeads.hg19)
library(RnBeads.hg38)
library(RnBeads)
library(GenomicRanges)
library(liftOver)
library(isva)
library(biomaRt)
library(RNOmni)
library(PCAtools)
library(tidyverse)


# summarise gene regions using cpgs, averages, and the pc method


datadir <- "/path/to/data"
savedir <- "/path/to/output"

dir.create(savedir, showWarnings=TRUE)

cleaned_str = "cleaned_global_"

get_pcs <- function(mvals){
    # estimate dims and get sig pcs for this gene
    head(mvals)
    dim(mvals)

    # normalize using inverse normal transform
    int_mvals <- apply(mvals, 1, RankNorm) %>% t()
    head(int_mvals)

    # get PCs
    pca_res <- pca(int_mvals)

    # get dim using gavish donoho
    # rows = variables
    # cols = observations
    gv_res <- chooseGavishDonoho(int_mvals, var.explained=pca_res$sdev^2, noise=1)
    gv_res

    # get dim using marchenko pastur
    mp_res <- chooseMarchenkoPastur(int_mvals, var.explained=pca_res$sdev^2, noise=1)
    mp_res

    print(paste("gavish donoho dim:", gv_res[[1]], "marchenko pastur dim:", mp_res[[1]]))

    names(pca_res)
    pcs <- pca_res$rotated
    head(pcs)

    est_dim <- gv_res[[1]]
    est_dim

    
    # subset pcs to estimated dimension
    sig_pcs <- pcs[,1:est_dim,drop=FALSE] %>%
        as.data.frame()
    #rownames(pcs) <- rownames(gene_methylation)
    head(sig_pcs)
    dim(sig_pcs)

    # get pct var explained
    eig_sq <- pca_res$sdev ^ 2
    pct_var <- eig_sq / sum(eig_sq)
    subset_pct_var <- pct_var[1:est_dim]

    # add the percent variance explained to the dataframe
    formatted_pcs <- rbind(subset_pct_var, sig_pcs)
    rownames(formatted_pcs)[1] <- "percent_variance"
    head(formatted_pcs)

    mp_var <- sum(pct_var[1:mp_res[[1]]])
    print(paste("Variance explained, GV:", round(sum(subset_pct_var), 2), "MP:", round(sum(mp_var), 2)))

    formatted_pcs
}

match_mvals <- function(mvals){
    # match the subject to those in the wgs filtered data
    head(mvals)
    dim(mvals)

    # subset to samples that we keep
    datafile <- paste0("../../output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds")
    keep_samples <- readRDS(datafile)
    head(keep_samples)
    dim(keep_samples)

    length(intersect(colnames(mvals), keep_samples$individualID))

    sub_mvals <- mvals[,keep_samples$individualID]
    head(sub_mvals)
    dim(sub_mvals)
    sum(is.na(sub_mvals))

    sub_mvals
}

main <- function(array_idx, region_idx, celltype_idx){
    # cell types to process
    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')

    # process each cell type seprately
    process_celltype <- function(cell_type){ 
        # load the normalized methylation data
        print(paste("Loading data for", cell_type, Sys.time()))
        datafile <- paste0(datadir, cleaned_str,
                           cell_type, "_lifted_hg38_mvals.rds")
        datafile
        mvals <- readRDS(datafile)
        head(mvals)
        dim(mvals)

        # subset mvals to samples that we keep
        matched_mvals <- match_mvals(mvals)

        # get the global pcs
        print(paste("Computing PCs for", cell_type, Sys.time()))
        pcs <- get_pcs(matched_mvals)
        head(pcs)

        # save global pc
        savefile <- paste0(savedir, cleaned_str, cell_type, "_global_meth_pcs.rds")
        savefile
        saveRDS(pcs, savefile)

        1
    }

    lapply(cell_types, process_celltype)
}
