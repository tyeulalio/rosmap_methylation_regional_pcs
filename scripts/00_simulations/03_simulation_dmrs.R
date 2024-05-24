library(RNOmni)
library(DMRcate)
library(minfi)
library(limma)
library(GenomicRanges)
library(tidyverse)

# run DMR analysis on simulated regions


savedir <- paste0("../../output/03_simulation_dmrs/")
dir.create(savedir, showWarnings = FALSE)

load_simulation <- function(num_sites, num_samples, 
                            percent_meth_difference, percent_sites_dm,
                            dmr_length, N
                            ){
    # load methylation for this set of parameters
    datadir <- paste0("../../output/00_simulated_data/")
    (filename <- paste0("simulated_data",
                        "_numSites", num_sites,
                        "_numSamples", num_samples,
                        "_pct_meth_diff", percent_meth_difference,
                        "_dmr_length", dmr_length,
                        "_pct_sites_dm", percent_sites_dm,
                        "_N", N))
    (datafile <- paste0(datadir, filename, ".rds"))

    # load the simulation data
    sims <- readRDS(datafile)
    names(sims)

    # select the raw data
    raw <- sims[['raw_data']]
    head(raw[[1]])
    length(raw)
    head(raw[[1]][[1]])

    # add the gene index number to each simulated region
    formatted_raw <- raw
    i <- 1000
    raw[[i]]
    for (i in seq_along(raw)){
        if (all(is.na(raw[[i]]))){
            formatted_raw[[i]] <- NA
            next
        }

        df1 <- raw[[i]][[1, drop=F]] %>%
            as.data.frame()
        df2 <- raw[[i]][[2, drop=F]] %>%
            as.data.frame()
        df1
        df2

        if (all(is.na(df1)) | all(is.na(df2))){
            formatted_raw[[i]] <- NA
        }

        # case where there is only 1 df2 row
        if (ncol(df2) == 1){
            df2 <- df2 %>%
                t() %>%
                as.data.frame()
            names(df2) <- c('dmr_row_start', 'dmr_chrpos', 'dmr_number_sites')
        }
        df2

        formatted_raw[[i]][[1]] <- df1 %>%
            mutate(gene_id = i)
        formatted_raw[[i]][[2]] <- df2 %>%
            mutate(gene_id = i)

        #if (i %% 100 == 0) print(paste("Processed", i, "genes out of 1,000"))
        #print(i)
    }
    head(formatted_raw[[1]])

    formatted_raw <- formatted_raw[!is.na(formatted_raw)]

    # get the combined data
    meth_df <- do.call(rbind, lapply(formatted_raw, function(x) x[[1]]))
    head(meth_df)[,1:10]
    dim(meth_df)

    dmr_info_df <- do.call(rbind, lapply(formatted_raw, function(x) x[[2]]))
    head(dmr_info_df)

    simulated_data <- list(meth_df=meth_df, dmr_info_df=dmr_info_df)
    simulated_data
}

format_simulation_data <- function(simulated_data, runname){
    # format the simulated data and take what we need
    names(simulated_data)

    # separate the coverage data from the meth percent data
    meth_df <- simulated_data[['meth_df']]
    head(meth_df)

    coverage_df <- meth_df %>%
        select(gene_id, chrpos:mprob.diff, 
               starts_with("SC"))
    head(coverage_df)

    percent_meth_df <- meth_df %>%
        select(gene_id, chrpos:mprob.diff, 
               starts_with("PM"))
    head(percent_meth_df)

    # save these dataframes
    formatted_sims <- list(coverage_df=coverage_df, 
                           percent_meth_df=percent_meth_df)
    (savefile <- paste0(savedir, 
                      "formatted_sims",
                      runname,
                      ".rds"))
    saveRDS(formatted_sims, savefile)

    formatted_sims
}

process_methylation <- function(long_meth){
    # process and normalize the methylation data
    # - remove low variance cpgs
    # - rank normalize
    head(long_meth)
    
    # make a cpg x sample data frame with meth values
    meth_df <- long_meth %>%
        select(chrpos, sample, PM, control) %>%
        mutate(control = ifelse(control, "control", "case")) %>%
        unite(sample, control, sample) %>%
        spread(sample, PM) %>%
        column_to_rownames('chrpos')
    head(meth_df)
    
    # remove zero variance cpgs
    # this is a standard practice
    var_meth <- meth_df[apply(meth_df, 1, var, na.rm=TRUE) != 0,]
    dim(var_meth)
    head(var_meth)
    nrow(meth_df) - nrow(var_meth)

    # apply inverse normal tranformation
    int_meth <- apply(var_meth, 1, RankNorm) %>%
        t() %>%
        as.data.frame()
    head(int_meth)
    
    return(int_meth)
}

normalize_methylation <- function(percent_meth_df, num_samples){
    # perform normalization
    # get long form of meth data
    long_meth <- percent_meth_df %>%
        gather('type', 'value', -gene_id, -chrpos, -dmsites, -mprob, -mprob.diff) %>%
        separate(type, c('type', 'sample'), sep="\\.") %>%
        spread(type, value) %>%
        # first half of samples are controls
        mutate(control = sample <= num_samples/2)
    head(long_meth)

    # make a cpg x sample data frame with meth values
    meth_df <- long_meth %>%
        select(gene_id, chrpos, sample, PM, control) %>%
        mutate(control = ifelse(control, "control", "case")) %>%
        unite(sample, control, sample) %>%
        spread(sample, PM) %>%
        mutate(gene_id = paste0("gene", gene_id)) %>%
        unite(pos, gene_id, chrpos, sep="_") %>%
        column_to_rownames('pos')
    head(meth_df)

    # remove zero variance cpgs
    # this is a standard practice
    var_meth <- meth_df[apply(meth_df, 1, var, na.rm=TRUE) != 0,]
    dim(var_meth)
    head(var_meth)
    nrow(meth_df) - nrow(var_meth)

    # apply inverse normal tranformation
    int_meth <- apply(var_meth, 1, RankNorm) %>%
        t() %>%
        as.data.frame()
    head(int_meth)

    # remove any gene with only a single CpG
    # these can't run through DMR
    genes <- data.frame(gene_pos = rownames(int_meth)) %>%
        separate(gene_pos, c('gene_id', 'pos'), sep='_', remove=FALSE)
    head(genes)

    single_genes <- genes %>%
        dplyr::count(gene_id) %>%
        filter(n == 1) %>%
        left_join(genes)
    single_genes

    keep_meth <- int_meth[!rownames(int_meth) %in% single_genes$gene_pos,]
    head(keep_meth)

    keep_meth
}

run_limma_dmp <- function(normalized_meth, pheno, runname, limma_savefile){
    # run limma dm analysis


    # set up the model + design matrix
    model_content <- paste0("~ status")
    design <- model.matrix(as.formula(model_content), data=pheno)
    
    # run linear model
    lmfit <- lmFit(normalized_meth, design)
    
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

    print("Saving limma results")
    saveRDS(results, limma_savefile)

    results
}


run_dmrcate <- function(formatted_sims, runname, num_samples, dmrcate_savefile, limma_savefile){

    # run dmrcate for dmr analysis
    names(formatted_sims)
    percent_meth_df <- formatted_sims[['percent_meth_df']]
    head(percent_meth_df)

    # normalize methylation
    normalized_meth <- normalize_methylation(percent_meth_df, num_samples)
    head(normalized_meth)

    # create phenotype dataframe
    pheno <- data.frame(samples=colnames(normalized_meth)) %>%
        mutate(status=str_remove(samples, "_[0-9]*")) %>%
        mutate(status=as.numeric(as.factor(status)))
    head(pheno)

    # run cpg-level DM analysis with limma
    limma_res <- run_limma_dmp(normalized_meth, pheno, runname, limma_savefile)
    head(limma_res)

    # run dmrcate 
    # select significance threshold
    fdr_thresh <- 0.05

    # create grset to make input for dmrcate
    head(normalized_meth)
    head(limma_res)

    meth_annotated <- limma_res %>%
        mutate(gene_pos = rownames(limma_res)) %>%
        separate(gene_pos, c('chr', 'pos'), sep='_') %>%
        select(stat=t, chr, pos, diff=logFC, fdr=bh_pval) %>%
        mutate(pos = as.numeric(pos))
    head(meth_annotated)

    (nsig <- sum(meth_annotated$fdr < fdr_thresh))

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
    str(dmr_results)

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

    print("Saving dmrcate results")
    saveRDS(dmr_df, dmrcate_savefile)

    1
}

main <- function(runnum){
    # run DMR analysis on the simulated regions
    num_sites_range <- c(20,50)
    #num_samples_range <- c(50,500,5000)
    num_samples_range <- c(50,500)
    #percent_meth_difference_range <- seq(0.00001, 0.1, 0.01)
    percent_meth_difference_range <- c(0)
    percent_sites_dm_range <- seq(0.25,1,0.25)

    # differential methylation parameters
    runs <- expand.grid(num_sites=num_sites_range,
                        num_samples=num_samples_range,
                        percent_meth_difference=percent_meth_difference_range,
                        percent_sites_dm=percent_sites_dm_range
    ) %>%
        mutate(run=row_number())
    dim(runs)
    head(runs)

    # for testing
    num_sites <- num_sites_range[2]
    num_samples <- num_samples_range[1]
    percent_meth_difference <- percent_meth_difference_range[1]
    percent_sites_dm <- percent_sites_dm_range[1]

    run <- runs[7,]
    process_run <- function(run){
        num_sites <- run[['num_sites']] %>% as.numeric()
        num_samples <- run[['num_samples']] %>% as.numeric()
        percent_meth_difference <- run[['percent_meth_difference']] %>% as.numeric()
        percent_sites_dm <- run[['percent_sites_dm']] %>% as.numeric()

        dmr_length <- num_sites / 2
        N <- 1000

        # process the run
        print(paste("Parameters:", "numsite:", num_sites, 
                    "numsamples:", num_samples,
                    "percent_meth_difference:", percent_meth_difference,
                    "percent_sites_dm:", percent_sites_dm,
                    "dmrlength:", dmr_length,
                    "N:", N,
                    "Run", run[['run']], "out of", nrow(runs)))

        (runname <- paste0(
                            "_numSites", num_sites,
                            "_numSamples", num_samples,
                            "_pct_meth_diff", percent_meth_difference,
                            "_dmr_length", dmr_length,
                            "_pct_sites_dm", percent_sites_dm,
                            "_N", N))

        # skip if output exists
        (dmrcate_savefile <- paste0(savedir,
                           "dmrcate_results",
                           runname,
                           ".rds"
                        ))
        if (file.exists(dmrcate_savefile)){
            print("DMRcate output exists. Skipping run.")
            return(1)
        }
        (limma_savefile <- paste0(savedir, 
                            "limma_dmp_results",
                            runname,
                            ".rds"
                            ))
        if (file.exists(limma_savefile)){
            print("limma output exists but DMRcate does not. Skipping run.")
            return(1)
        }

        # need to perform CpG-level DM first
        # create a dataframe of all simulated sites
        message("Loading simulated methylation")
        simulated_data <- load_simulation(num_sites, num_samples, percent_meth_difference, percent_sites_dm, dmr_length, N)

        # format simulation data
        message("Formatting simulated data")
        formatted_sims <- format_simulation_data(simulated_data, runname)

        # try DMRcate analysis 
        message("Running limma and DMRcate")
        run_dmrcate(formatted_sims, runname, num_samples, dmrcate_savefile, limma_savefile)

        1
    }

    #process_run(runs[runnum,])

    res <- apply(runs, 1, process_run)
}

# get command line arguments and run main
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
# run number
runnum = args[1]

main(runnum)
