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
# setwd("D:/montgomery_lab/deconvolution_project/rosmap_methylation_regional_pcs/scripts/00_simulations/")
source("RRBSsimulation.R")
source("RRBSmodeling.R")

# create a directory to save the data
savedir <- paste0("../../output/00_simulated_data/")
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
    
    summarized_region <- list(avgs=avgs, rpcs=rpcs, full_rpcs=pc_res)
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
                                     res <- tryCatch({get_region(region_name,
                                                num_sites, 
                                                num_samples, 
                                                dmr_length,
                                                percent_meth_difference, 
                                                percent_sites_dm)
                                            },
                                            error = function(e){
                                                print(paste("Error:", e))
                                                return(list(avgs=NA, rpcs=NA, raw_data=NA))
                                            }
                                     )
                                }
                            )

    if (all(is.na(summarized_regions))){
        return(NA)
    }
    print(names(summarized_regions[[1]]))

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
    
    rpcs_full <- lapply(summarized_regions, function(x) x$full_rpcs)
    head(rpcs_full[[1]])
    
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
                          rpcs_full=rpcs_full,
                          parameters=parameters)
    simulated_dat
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

main <- function(runnum, N){
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
    # num_sites = 50
    # num_samples = 500
    # percent_meth_difference = 0.1
    # dmr_length = 25
    # percent_sites_dm = 0.1

    
    # number of regions simulate
    # N = 5
    
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
    }
        
        
    # parameter ranges that are being tested
    #num_sites_range <- c(20,50)
    #num_samples_range <- c(50,500,5000)
    #percent_sites_dm_range <- seq(0,0.75,0.25)
    #percent_meth_difference_range <- seq(0.1, 0.9, 0.1)

    # paramters to test on second run
    num_sites_range <- c(20,50)
    num_samples_range <- c(50,500)
    percent_sites_dm_range <- seq(0,0.75,0.25)
    percent_meth_difference_range <- seq(0.01, 0.2, 0.01)
    dmr_length_range = c(250, 2000)

    runs <- expand.grid(num_sites=num_sites_range,
                        num_samples=num_samples_range,
                        percent_meth_difference=percent_meth_difference_range,
                        percent_sites_dm=percent_sites_dm_range
    ) %>%
        mutate(run=row_number())
    dim(runs)
    head(runs)

    
    process_run <- function(run){
        num_sites <- run[['num_sites']] %>% as.numeric()
        num_samples <- run[['num_samples']] %>% as.numeric()
        percent_meth_difference <- run[['percent_meth_difference']] %>% as.numeric()
        percent_sites_dm <- run[['percent_sites_dm']] %>% as.numeric()

        #dmr_length <- num_sites / 10

        print(paste("Processing run", run[['run']], "out of", nrow(runs)))
        print(paste("Parameters:", "numsite:", num_sites, 
                    "numsamples:", num_samples,
                    "percent_meth_difference:", percent_meth_difference,
                    "percent_sites_dm:", percent_sites_dm,
                    "dmrlength:", dmr_length,
                    "N:", N))
        
        res <- run_simulation(num_sites, num_samples, percent_meth_difference, dmr_length, percent_sites_dm, N)
        return(1)
    }


    #res <- apply(runs, 1, process_run)
    process_run(runs[runnum,])
    return(1)
}


# get command line arguments and run main
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
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
