library(ggplot2)
library(tidyverse)

# example script that processes the simulation results from simulation_test.R script
# this script will look at the results from 1 parameter change in the simulations


# check one result file to explore the structure
check_example_result <- function(){
    num_sites <- 50 
    num_samples <- 500
    percent_meth_difference <- 0.5
    percent_sites_dm <- 0.5
    N <- 1000
    dmr_length <- num_sites / 2

    # name of output files are formatted like this
    datadir <- paste0("./simulated_data/")
    (filename <- paste0("simulated_data",
                        "_numSites", num_sites,
                        "_numSamples", num_samples,
                        "_pct_meth_diff", percent_meth_difference,
                        "_dmr_length", dmr_length,
                        "_pct_sites_dm", percent_sites_dm,
                        "_N", N))
    (datafile <- paste0(datadir, filename, ".rds"))

    simulation_run <- readRDS(datafile)
    names(simulation_run)

    test <- simulation_run$parameters
    class(test)
    names(test)
    head(test)
    length(test)
    
}

# main function to run the script
main <- function(){
   # meth simulation parameters
    # num_sites = number of CpG sites to simulate in the region
    # num_samples = number of samples to simulate data for
    # percent_meth_difference = how different should controls and cases be 
    # - at DMR sites
    # dmr_length = length of differentially methylated regions 
    # - this is not the same as # of differentially methylated CpGs
    # percent_sites_dm = percentages of CpGs that should be differentially
    # methylated 

    # ranges for parameters tested:
    # every combination of parameters was run
    #  num_sites_range <- c(20,50)
    #  num_samples_range <- c(50,500,5000)
    #  percent_meth_difference_range <- seq(0.1, 0.9, 0.1)
    #  percent_sites_dm_range <- seq(0,1,0.25)

    # parameters that do not change
    # N = 1000 for everything (number of regions simulated)
    # dmr_length = num_sites / 2

    check_example_result()
}

main()
