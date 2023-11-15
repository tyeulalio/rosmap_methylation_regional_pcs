## Adapted some of this code from the example provided by the simulation model.
## This code runs the simulation analysis for regional pcs.
## See main function at the bottom for the main workflow.

library(regionalpcs)
library(RNOmni)
library(ggplot2)
library(limma)
library(tidyverse)

# *** Set the working drive to wherever you have the modeling scipts ***
## First source both the "RRBSsimulation.R" and "RRBSmodeling.R" scripts and make sure that all required libraries are installed
## Install the appropriate packages, some are out of date, you may get
## warning messages about invalid components but the code will still run
#setwd("D:/OneDrive/MontgomeryLab/deconvolution/RRBS-sim/")
source("RRBSsimulation.R")
source("RRBSmodeling.R")

# create a directory to save the data
savedir <- paste0("simulated_data/")
dir.create(savedir)



simulate_region <- function(num_sites, num_samples, dmr_length,
                          percent_meth_difference, percent_sites_dm){
    # simulate a single gene's methylation
    region_sim <- sim.bed.vg(nsites=num_sites,
               nsamp=num_samples,
               dm.len=dmr_length,
               prop.diff=percent_meth_difference,
               pm.diff=percent_sites_dm)
    region_sim
    
    return(region_sim)
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


summarize_region <- function(processed_meth, region_name){
    # summarize region using avgs and regional pcs
    head(processed_meth)
    
    # compute averages
    avgs <- colMeans(processed_meth) %>%
        t()
    rownames(avgs) <- region_name
    head(avgs)
    
    # -- compute regional pcs
    
    # create a region map for the function
    region_map <- data.frame(region=region_name, cpg=rownames(processed_meth))
    head(region_map)
    
    # compute regional pcs using the regionalpcs package
    pc_res <- compute_regional_pcs(processed_meth, region_map, pc_method='gd')
    names(pc_res)
    
    # probably save the pcs and averages at this step
    # may want to look back at percent variance/loadings later
    
    rpcs <- pc_res$regional_pcs
    head(rpcs)
    dim(rpcs)
    
    summarized_region <- list(avgs=avgs, rpcs=rpcs)
    summarized_region
}

get_region <- function(region_name, num_sites, num_samples, dmr_length,
                     percent_meth_difference, percent_sites_dm){

    
    # simulate a single region's methylation
    region_sim <- simulate_region(num_sites, num_samples, dmr_length,
                                    percent_meth_difference, percent_sites_dm)
    
    # may want to save region_sim at this point to look back at later
    
    ##View the first few rows of simulated data.  
    # "dmsites" in an indicator taking values (-1,0,1), 
    # "mprob" is the PM for group 1, "mprob.diff" is the PM for group 2.  
    # The next 8 columns contain the scores and 
    # the final 8 contain the PM values, with 1-4 corresponding to group 1 and 5-8
    # corresponding to group 2. 
    head(region_sim[[1]])
    
    simulated_meth <- region_sim[[1]]
    dim(simulated_meth)
    head(simulated_meth)
    
    # tells you direction of DM difference for controls vs cases
    # table(simulated_meth$dmsites)
    
    
    # get long form of meth data
    long_meth <- simulated_meth %>%
        gather('type', 'value', -chrpos, -dmsites, -mprob, -mprob.diff) %>%
        separate(type, c('type', 'sample'), sep="\\.") %>%
        spread(type, value) %>%
        # first half of samples are controls
        mutate(control = sample <= num_samples/2)
    head(long_meth)
    
    
    ##View the DMRs.  
    # Column 1 contains the row indices for the DMR starts, 
    # column 2 contains the chrpos values, and 
    # column 3 contains the number of sites
    head(region_sim[[2]])
    if (!is.null(nrow(region_sim[[2]]))){
        colnames(region_sim[[2]]) <- c('dmr_row_start', 'dmr_chrpos', 'dmr_number_sites')
    }
       
    
    # process and normalize the methylation data
    processed_meth <- process_methylation(long_meth)
    
    # summarize region using avgs and regional pcs
    summarized_region <- summarize_region(processed_meth, region_name)
    names(summarized_region)
    
    summarized_region$raw_data <- region_sim
    
    summarized_region
}

get_N_regions <- function(num_sites, num_samples, percent_meth_difference,
                          dmr_length, percent_sites_dm, N){
    # get N regions with specific parameters for simulation meth
    # summarized as averages and regional pcs
   
    # meth simulation parameters
    # num_sites = number of CpG sites to simulate in the region
    # num_samples = number of samples to simulate data for
    # percent_meth_difference = how different should controls and cases be 
    # - at DMR sites
    # dmr_length = length of differentially methylated regions 
    # - this is not the same as # of differentially methylated CpGs
    # percent_sites_dm = percentages of CpGs that should be differentially
    # methylated
    # num_sites = 500
    # num_samples = 500
    # percent_meth_difference = 0.5
    # dmr_length = 250
    # percent_sites_dm = 0.1
    
    # number of regions simulate
    # N = 1000
    
    # summarize N regions
    summarized_regions <- lapply(seq_len(N), 
                                 function(x) {
                                     message(paste("Generating region", x))
                                     region_name <- paste0("region", x)
                                     get_region(region_name,
                                                num_sites, 
                                                num_samples, 
                                                dmr_length,
                                                percent_meth_difference, 
                                                percent_sites_dm)
                                     }
                                 )
    summarized_regions[[1]]
    
    # create summary type data frames
    avgs <- lapply(summarized_regions, function(x) x$avgs) %>%
        do.call(rbind, .)
    head(avgs)
    
    rpcs <- lapply(summarized_regions, function(x) x$rpcs) %>%
        do.call(rbind, .)
    head(rpcs)
    
    # store the raw simulated data also
    head(summarized_regions[[1]]$raw_data)
    raw_data <- lapply(summarized_regions, function(x) x$raw_data)
    head(raw_data[[1]])
    
    # store the parameters
    parameters <- list(num_sites=num_sites, 
                       num_samples=num_samples,
                       percent_meth_different=percent_meth_difference,
                       dmr_length=dmr_length,
                       percent_sites_dm=percent_sites_dm,
                       N=N)
    
    # combine the results and parameters together
    simulated_dat <- list(avgs=avgs, rpcs=rpcs,
                          raw_data=raw_data,
                          parameters=parameters)
    simulated_dat
}

#meth <- avgs # for testing
run_dmp <- function(meth, pheno){
    # run the dmp analysis
    head(meth) 
    head(pheno)
    
    stopifnot(colnames(meth) == pheno$samples)
    
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
    
    return(results)
}

## -- checking on simulated data
check_data <- function(){
    simulated_meth <- simulated_dat$raw_data[[1]][[1]]
    dim(simulated_meth)
    head(simulated_meth)[,1:10]
    
    dmr_locations <- simulated_dat$raw_data[[1]][[2]]
    head(dmr_locations)
    
    # tells you direction of DM difference for controls vs cases
    # table(simulated_meth$dmsites)
    
    # get long form of meth data
    long_meth <- simulated_meth %>%
        gather('type', 'value', -chrpos, -dmsites, -mprob, -mprob.diff) %>%
        separate(type, c('type', 'sample'), sep="\\.") %>%
        spread(type, value) %>%
        # first half of samples are controls
        mutate(control = sample <= num_samples/2) 
    head(long_meth)
    
    # where are the DMRs?
    ggplot(simulated_meth) +
        geom_point(aes(x=chrpos, y=dmsites)) 
    
    # long_meth %>%
    #     group_by(chrpos, control, dmsites) %>%
    #     summarize(pm_mean = mean(PM),
    #               pm_median = median(PM)) %>%
    #     select(-pm_median) %>%
    #     spread(control, pm_mean) %>%
    #     mutate(diff = `FALSE` - `TRUE`) %>%
    #     group_by(dmsites) %>%
    #     summarize(mean=mean(diff),
    #               median=median(diff))
    
    # just visualizing simulated meth
    ggplot(long_meth, aes(x=as.factor(chrpos), y=PM, color=control)) +
        geom_boxplot() +
        facet_wrap(. ~ dmsites, ncol=1, scales='free_x') +
        ggtitle(filename)
    
    ggplot(long_meth) +
        geom_density(aes(x=PM, color=control)) +
        facet_wrap(. ~ dmsites, scales='free_y')
    
    ggplot(long_meth) +
        geom_histogram(aes(x=PM, fill=control)) +
        facet_wrap(control ~ dmsites, scales='free_y')
}

main <- function(){
    # simulate the data
    # meth simulation parameters
    # num_sites = number of CpG sites to simulate in the region
    # num_samples = number of samples to simulate data for
    # percent_meth_difference = how different should controls and cases be 
    # - at DMR sites
    # dmr_length = length of differentially methylated regions 
    # - this is not the same as # of differentially methylated CpGs
    # percent_sites_dm = percentages of CpGs that should be differentially
    # methylated
    num_sites = 50
    num_samples = 500
    percent_meth_difference = 0.1
    dmr_length = 25
    percent_sites_dm = 0.1
    
    # number of regions simulate
    N = 100
    
    # simulate and run DM analysis
    run_simulation <- function(num_sites, num_samples,
                               percent_meth_difference,
                               dmr_length, percent_sites_dm, N){
        
        # store the simulated data
        (filename <- paste0("simulated_data",
                            "_numSites", num_sites,
                            "_numSamples", num_samples,
                            "_pct_meth_diff", percent_meth_difference,
                            "_dmr_length", dmr_length,
                            "_pct_sites_dm", percent_sites_dm,
                            "_N", N))
        (savefile <- paste0(savedir, filename, ".rds"))
        
        # load the output if we ran these parameters before
        # otherwise, simulate the data
        if (file.exists(savefile)){
            simulated_dat <- readRDS(savefile)
        } else{
            # simulate the region's meth for N regions
            simulated_dat <- get_N_regions(num_sites, num_samples,
                                           percent_meth_difference,
                                           dmr_length, percent_sites_dm, N)
            
            # store the simulated data
            (filename <- paste0("simulated_data",
                                "_numSites", num_sites,
                                "_numSamples", num_samples,
                                "_pct_meth_diff", percent_meth_difference,
                                "_dmr_length", dmr_length,
                                "_pct_sites_dm", percent_sites_dm,
                                "_N", N))
            (savefile <- paste0(savedir, filename, ".rds"))
            saveRDS(simulated_dat, savefile)
        }
        
        
        # grab normalized methylation data
        rpcs <- simulated_dat$rpcs
        avgs <- simulated_dat$avgs
        
        print(paste("number of rpcs:", nrow(rpcs)))
        print(paste("number of avgs:", nrow(avgs)))
        
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
        avgs_sig <- get_sig_counts(avgs_res)
        rpcs_sig <- get_sig_counts(rpcs_res)
        
        tibble(num_sites=num_sites,
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
    }
    
    
    ## -- Example test for changing percent meth
    (test_range <- seq(0.1, 1, 0.1))
    res <- lapply(test_range, function(x) run_simulation(num_sites, num_samples,
                                                   x,
                                                   dmr_length, percent_sites_dm, N))
    res_df <- do.call(rbind, res)
    head(res_df)
    
    # just a quick view of the results
    # get long version
    long_res <- res_df %>%
        select(percent_meth_difference, bf_sig_avgs, bf_sig_rpcs) %>%
        gather('sig_type', 'num_sig', bf_sig_avgs, bf_sig_rpcs) %>%
        mutate(summary_type = str_remove(sig_type, "bf_sig_"))
    head(long_res)
    
    ggplot(long_res) +
        geom_point(aes(x=percent_meth_difference,
                       y=num_sig,
                       color=summary_type),
                   size=2) +
        theme_bw() +
        theme(text=element_text(size=20)) +
        xlab(paste("Percent methylation difference")) +
        ylab(paste("Number of significant DM results")) +
        guides(color=guide_legend(title="Summary type"))
        
}

main()











