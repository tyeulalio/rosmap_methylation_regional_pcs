# Run differential methylation analysis on the simulated regions

#library(regionalpcs)
library(RNOmni)
library(ggplot2)
library(limma)
library(tidyverse)


# create a directory to save the data
savedir <- paste0("../../output/01_dmp_results/")
dir.create(savedir)

datadir <- paste0("simulated_data/")



#meth <- avgs # for testing
run_dmp <- function(meth, pheno){
    # run the dmp analysis
    head(meth) 
    head(pheno)

    stopifnot(colnames(meth) == pheno$samples)

    # store results here
    res_list <- list()


    # store dimensions of the full dataframe
    res_list$full_meth_dims <- dim(na.omit(meth))

    # remove rows with NA values - failed simulations
    full_meth <- na.omit(meth)

    # apply inverse normal tranformation
    int_meth <- apply(full_meth, 1, RankNorm) %>%
        t() %>%
        as.data.frame()
    head(int_meth)
    
    meth <- int_meth
    
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

    res_list$results <- results
    
    return(res_list)
}

main <- function(runnum, N){
    
    # run DM analysis
    run <- runs[1,] # for testing interactively
    process_run <- function(run){
        num_sites <- run[['num_sites']] %>% as.numeric()
        num_samples <- run[['num_samples']] %>% as.numeric()
        percent_meth_difference <- run[['percent_meth_difference']] %>% as.numeric()
        percent_sites_dm <- run[['percent_sites_dm']] %>% as.numeric()

        dmr_length <- num_sites / 2

        print(paste("Processing run", run[['run']], "out of", nrow(runs)))
        print(paste("Parameters:", "numsite:", num_sites, 
                    "numsamples:", num_samples,
                    "percent_meth_difference:", percent_meth_difference,
                    "percent_sites_dm:", percent_sites_dm,
                    "dmrlength:", dmr_length,
                    "N:", N))
        
        # store the simulated data
        (filename <- paste0("simulated_data",
                            "_numSites", num_sites,
                            "_numSamples", num_samples,
                            "_pct_meth_diff", percent_meth_difference,
                            "_dmr_length", dmr_length,
                            "_pct_sites_dm", percent_sites_dm,
                            "_N", N))
        (datafile <- paste0(datadir, filename, ".rds"))
        
        # load the output if we ran these parameters before
        # otherwise, simulate the data
        if (file.exists(datafile)){
            results <- readRDS(datafile)
        } else{
            print("Results for this run do not exist")
            return(1)
        }
        
        head(results)

        # grab the averages and rpcs
        avgs <- results$avgs
        rpcs <- results$rpcs

        # check if we can detect differential methylation
        # maybe do this as a block of about 15k simulated regions
        # like standard DM analysis? Apply multiple test correction too.
        # make a phenotype dataframe
        pheno <- data.frame(samples=colnames(avgs)) %>%
            mutate(status=str_remove(samples, "_[0-9]*")) %>%
            mutate(status=as.numeric(as.factor(status)))
        head(pheno)
        
        avgs_res <- run_dmp(avgs, pheno)
        rpcs_res <- run_dmp(rpcs, pheno)


        dm_res <- list(avgs_res=avgs_res, rpcs_res=rpcs_res)
        
        # compare methods
        head(avgs_res)
        
        dm_res <- avgs_res
        sig_thresh <- 0.05
        get_sig_counts <- function(dm_res){
            bf_sig <- sum(dm_res$bf_pval < sig_thresh)
            bh_sig <- sum(dm_res$bh_pval < sig_thresh)
            uc_sig <- sum(dm_res$P.Value < sig_thresh)
            list(bf_sig=bf_sig,
                 bh_sig=bh_sig,
                 uc_sig=uc_sig)
        }
        avgs_sig <- get_sig_counts(avgs_res$results)
        rpcs_sig <- get_sig_counts(rpcs_res$results)
        
        num_sig <- tibble(num_sites=num_sites,
               num_samples=num_samples,
               percent_meth_difference=percent_meth_difference,
               dmr_length=dmr_length,
               percent_sites_dm=percent_sites_dm,
               N=N,
               num_rpcs=nrow(rpcs),
               
               bf_sig_avgs=avgs_sig$bf_sig,
               bf_sig_rpcs=rpcs_sig$bf_sig,
               
               bh_sig_avgs=avgs_sig$bh_sig,
               bh_sig_rpcs=rpcs_sig$bh_sig,
               
               uc_sig_avgs=avgs_sig$uc_sig,
               uc_sig_rpcs=rpcs_sig$uc_sig)
        dm_res$num_sig <- num_sig

        (filename <- paste0("DM_results",
                            "_numSites", num_sites,
                            "_numSamples", num_samples,
                            "_pct_meth_diff", percent_meth_difference,
                            "_dmr_length", dmr_length,
                            "_pct_sites_dm", percent_sites_dm,
                            "_N", N))
        (savefile <- paste0(savedir, filename, ".rds"))
        saveRDS(dm_res, savefile)
    }
    
    
    num_sites_range <- c(20,50)
    num_samples_range <- c(50,500,5000)
    percent_meth_difference_range <- seq(0.1, 0.9, 0.1)
    percent_sites_dm_range <- seq(0,1,0.25)

    runs <- expand.grid(num_sites=num_sites_range,
                        num_samples=num_samples_range,
                        percent_meth_difference=percent_meth_difference_range,
                        percent_sites_dm=percent_sites_dm_range
    ) %>%
        mutate(run=row_number())
    dim(runs)
    head(runs)
    

    # for testing
    #runnum <- 1
    #process_run(runs[runnum,])

    res <- apply(runs, 1, process_run)
    return(1)
}


# get command line arguments and run main
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).", call.=FALSE)
} 
# run number
runnum = args[1]

# number of regions simulate
N = 1000

start <- Sys.time()
main(runnum, N)
end <- Sys.time()

print(paste("start time:", start))
print(paste("end time:", end))
print(paste("duration:", end-start))
