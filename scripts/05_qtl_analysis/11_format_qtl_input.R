library(tidyverse)

summary_types <- c('avgs', 'pcs', 'cpgs')
cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')

chromdir = ""

datadir <- "/path/to/data/"
savedir <- "/path/to/output/"
without_proportions=FALSE
cleaned_str = "cleaned_global_"

dir.create(savedir, showWarnings=FALSE)


match_covariates <- function(cell_type, region_type){
    # match up the covariates for the qtl analysis
    # ancestry - genotype PCs
    # global methylation PCs
    # other phenotype variables
    # cell type composition (for bulk only)

    orig_dir <- "/path/to/original/genotype/data/"

    # load genotype PCs
    (datafile <- paste0(orig_dir, "output/05_qtl_analysis/08_format_eigenstrat/", chromdir, "genotype_pcs_hg38.rds"))
    geno_pcs <- readRDS(datafile)
    head(geno_pcs)

    # load global meth pcs
    datafile <- paste0("../../output/05_qtl_analysis/10_global_pcs/", cleaned_str, cell_type, "_global_meth_pcs.rds")
    meth_pcs <- readRDS(datafile)
    head(meth_pcs)

    # load phenotype data
    datafile <- paste0(datadir, 
                       cell_type, "_formatted_clinical.rds")
    clinical <- readRDS(datafile)
    head(clinical)

    # combine covariates together
    # xqtl paper corrected for: RIN, PMI, study index, age death, sex
    # make all covars numeric
    head(clinical)
    covars <- c('pmi', 'Study', 'age_death', 'msex', 'meth.batch', 'neuron', 'astro', 'endo', 'oligo_opc')
    if (without_proportions){
        covars <- c('pmi', 'Study', 'age_death', 'msex', 'meth.batch')
    }
    covars
    sub_clinical <- clinical %>%
        select(wgs.specimenID, individualID, all_of(covars)) %>%
        mutate(studyROS = as.numeric(Study == 'ROS')) %>%
        select(-Study) 
    head(sub_clinical)

    head(geno_pcs)
    head(meth_pcs)

    # change names of meth pcs and genotype pcs
    colnames(meth_pcs) <- paste0("methylation_", colnames(meth_pcs))
    colnames(geno_pcs)[-1] <- paste0("genotype_", colnames(geno_pcs)[-1])


    # combine meth pcs, covariates, and geno pcs
    formatted_covariates <- meth_pcs %>%
        rownames_to_column('individualID') %>%
        filter(individualID != 'percent_variance') %>%
        inner_join(sub_clinical) %>%
        rename(specimenID=wgs.specimenID) %>%
        #select(-meth_specimenID) %>%
        inner_join(geno_pcs) %>%
        select(individualID, everything())
    head(formatted_covariates)

    dim(formatted_covariates)
    dim(geno_pcs)

    # fill in NA values
    full_cov <- formatted_covariates %>%
        mutate(pmi = ifelse(is.na(pmi), mean(pmi, na.rm=TRUE), pmi)) %>%
        mutate(msex = as.numeric(msex),
               meth.batch = as.numeric(meth.batch)
        )
    head(full_cov)
    sum(is.na(full_cov))

    full_cov
}

match_covariates_outcome <- function(covariates, cell_type, region_type, summary_type){
    # match the covariates with the outcome variable
    # outcome traits = summarized regional methylation values
    head(covariates)

    # load outcome 
    if (summary_type == 'cpgs'){ 
        datafile <- paste0("../../output/05_qtl_analysis/03_format_pc_beds/", cleaned_str, cell_type,
                           "_", "genome", "_formatted_", summary_type, ".bed.gz")
    } else{
        datafile <- paste0("../../output/05_qtl_analysis/03_format_pc_beds/", cleaned_str, cell_type,
                           "_", region_type, "_formatted_", summary_type, ".bed.gz")
    }
    datafile
    outcome_df <- read_tsv(datafile, progress=FALSE)
    head(outcome_df)
    dim(outcome_df)

    # hand the different summary types
    if (summary_type == 'cpgs'){
        head(covariates)
        head(outcome_df)

        common_individuals <- intersect(covariates$individualID, colnames(outcome_df)[5:ncol(outcome_df)])
        length(common_individuals)

        # match the two dataframes
        matched_covariates <- covariates[match(common_individuals, covariates$individualID),]
        head(matched_covariates)
        dim(matched_covariates)

        head(outcome_df)
        matched_outcome <- outcome_df %>%
            select('#chr', start, end, gene_id, all_of(common_individuals))
        head(matched_outcome)
        dim(matched_outcome)

        # change the colnames to specimenID to match geno data
        identical(matched_covariates$individualID, colnames(matched_outcome)[5:ncol(matched_outcome)])

        colnames(matched_outcome)[5:ncol(matched_outcome)] <- matched_covariates$specimenID
        head(matched_outcome)

        identical(matched_covariates$specimenID, colnames(matched_outcome)[5:ncol(matched_outcome)])

        # format the matched covariates
        matched_covariates <- matched_covariates %>%
            select(-individualID) %>%
            select(specimenID, everything())
        head(matched_covariates)
            
    }else{
        # regionalpcs and averages

        # get common samples
        common_samples <- intersect(covariates$specimenID, colnames(outcome_df)[5:ncol(outcome_df)])
        length(common_samples)

        # match the two dataframes
        matched_covariates <- covariates[match(common_samples, covariates$specimenID),] %>%
            select(-individualID) %>%
            select(specimenID, everything())
        head(matched_covariates)
        dim(matched_covariates)

        matched_outcome <- outcome_df %>%
            select('#chr', start, end, gene_id, all_of(common_samples))
        head(matched_outcome)[,1:10]
        dim(matched_outcome)
    }

    # make sure these match
    identical(matched_covariates$specimenID, colnames(matched_outcome)[5:ncol(matched_outcome)])

    # save the files
    head(matched_covariates)
    savefile <- paste(cell_type, region_type, summary_type, "matched_covariates.tsv", sep='_')
    savepath <- paste0(savedir, savefile)
    savepath
    write_tsv(matched_covariates, savepath, progress=FALSE)

    head(matched_outcome)
    savefile <- paste(cell_type, region_type, summary_type, "matched_pheno.bed", sep='_')
    savepath <- paste0(savedir, savefile)
    savepath
    write_tsv(matched_outcome, savepath, progress=FALSE)

    # gzip the output
    if (file.exists(paste0(savepath, ".gz"))) system(paste0("rm ", savepath, ".gz"))
    system(paste("gzip", savepath))

    1
}

main <- function(){
    # define the region types to test
    region_types <- c('full_gene')

    # get every combination of parameters to run
    runs <- expand.grid(summary_type = summary_types, 
                        cell_type=cell_types, 
                        region_type=region_types) %>%
        mutate(row_number = row_number())
    runs

    # run each combination of parameters
    run_all <- function(row){
        cell_type <- row[['cell_type']]
        summary_type <- row[['summary_type']]
        region_type <- row[['region_type']]
        rn <- row[['row_number']]

        print(paste("Processing", cell_type, region_type, summary_type, rn))

        # output file
        savefile <- paste(cell_type, region_type, summary_type, "matched_pheno.bed.gz", sep='_')
        savepath <- paste0(savedir, cleaned_str, savefile)
        savepath

        # skip if we processed this run already
        if (file.exists(savepath)){
            print(paste("OUTPUT FILE EXISTS. SKIPPING."))
            return(1)
        }

        # match phenotype and covariate files
        covariates <- match_covariates(cell_type, region_type)
        head(covariates)

        # match covariates with outcome trait 
        match_covariates_outcome(covariates, cell_type, region_type, summary_type)
        1
    }

    head(runs)
    apply(runs, 1, run_all)

}

main()
