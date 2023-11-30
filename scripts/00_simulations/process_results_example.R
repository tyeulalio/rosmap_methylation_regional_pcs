library(ggplot2)
library(tidyverse)

# example script that processes the simulation results from simulation_test.R script
# this script will look at the results from 1 parameter change in the simulations

savedir <- paste0("../../output/process_results_example/")
dir.create(savedir)

# use the colors defined in the color panels R script
source("../color_panels.R")

# checking how meth difference and percent site dm affects results
check_meth_difference <- function(num_sites, num_samples){

    # will vary percent_meth_difference and percent_sites_dm
    percent_meth_difference_range <- seq(0.1, 0.9, 0.1)
    percent_sites_dm_range <- seq(0,0.75,0.25)

    runs <- expand.grid(percent_meth_difference=percent_meth_difference_range,
                        percent_sites_dm=percent_sites_dm_range
    ) %>%
        mutate(rownum=row_number()) 
    head(runs)

    # fix these for this analysis
    #num_sites <- 50 
    #num_samples <- 500
    N <- 1000
    dmr_length <- num_sites / 2

    row <- runs[1,] # for testing
    load_data <- function(row){
        percent_meth_difference <- row[['percent_meth_difference']] %>% as.numeric()
        percent_sites_dm <- row[['percent_sites_dm']] %>% as.numeric()

        # name of output files are formatted like this
        datadir <- paste0("../../output/01_dmp_results/")
        (filename <- paste0("DM_results",
                            "_numSites", num_sites,
                            "_numSamples", num_samples,
                            "_pct_meth_diff", percent_meth_difference,
                            "_dmr_length", dmr_length,
                            "_pct_sites_dm", percent_sites_dm,
                            "_N", N))
        (datafile <- paste0(datadir, filename, ".rds"))

        # add this check for missing files
        # some runs may have failed - will look into these later
        if (!file.exists(datafile)){
            print("datafile does not exist")
            return(NA)
        }

        results <- readRDS(datafile)
        names(results)

        # grab the summarized results 
        summarized_res <- results$num_sig
        head(summarized_res)

        # add on the dimensions of the data frame
        summarized_res['rpcs_nregions'] <- results$rpcs_res$full_meth_dims[[1]]
        summarized_res['rpcs_nsamples'] <- results$rpcs_res$full_meth_dims[[2]]
        summarized_res['avgs_nregions'] <- results$avgs_res$full_meth_dims[[1]]
        summarized_res['avgs_nsamples'] <- results$avgs_res$full_meth_dims[[2]]

        summarized_res
    }

    results_list <- apply(runs, 1, load_data)
    results_df <- do.call(rbind, results_list) %>%
        na.omit()
    head(results_df)
    dim(results_df)

    # checking how many runs failed, not too many, 16 at most
    summary(results_df$rpcs_nregions)
    summary(results_df$avgs_nregions)

    # compute a proportion of significant results
    results_df <- results_df %>%
        mutate(bh_sig_prop_avgs = bh_sig_avgs / avgs_nregions,
               bh_sig_prop_rpcs = bh_sig_rpcs / rpcs_nregions,
               bf_sig_prop_avgs = bf_sig_avgs / avgs_nregions,
               bf_sig_prop_rpcs = bf_sig_rpcs / rpcs_nregions
        )
    head(results_df)

    # get a long version of data frame
    long_results <- results_df %>%
        select(percent_sites_dm,
               percent_meth_difference,
               bh_sig_prop_avgs:bf_sig_prop_rpcs
        ) %>%
        gather('prop_type', 'prop', bh_sig_prop_avgs:bf_sig_prop_rpcs) %>%
        mutate(prop_type = str_remove(prop_type, '_sig_prop')) %>%
        separate(prop_type, c('adjustment', 'summary_type'), sep='_') %>%
        filter(adjustment == 'bf')
    head(long_results)

    #num_sites <- 50 
    #num_samples <- 500
    #N <- 1000
    #dmr_length <- num_sites / 2
    (filename <- paste0(
                        "numSites", num_sites,
                        "_numSamples", num_samples,
                        #"_pct_meth_diff", percent_meth_difference,
                        "_dmr_length", dmr_length,
                        #"_pct_sites_dm", percent_sites_dm,
                        "_N", N))

    # plot the results
    p <- long_results %>%
        # format summary type to match colors map
        mutate(formatted_st = ifelse(summary_type == 'avgs', 'Avgs', 'rPCs')) %>%
        ggplot( 
                aes(x=percent_meth_difference,
                    y=prop,
                    color=formatted_st) ) +
        geom_point() +
        geom_line() +
        facet_wrap(. ~ percent_sites_dm) +
        scale_x_continuous(breaks=seq(0.1, 0.9, 0.1)) +
        scale_color_manual(values=summary_type_colors, name="Summary type") + # comes from the color panel file we sourced
        theme_bw() +
        theme(
              text=element_text(size=14),
              axis.text.x=element_text(angle=30, hjust=1, vjust=1),
              # remove minor grid lines
              panel.grid.minor=element_blank(),
              #panel.grid.major=element_blank(), # will do this later for paper
              strip.background=element_blank(), # remove facet strip background
              panel.border=element_rect(color='black',
                                        fill=NA,
                                        linewidth=1
                                        )
              ) +
        ylab("Proportion of regions with significant DM detected") +
        ggtitle(str_replace_all(filename, "_", " "))
    (savefile <- paste0(savedir, 
                        "percentSitesDm_percentMethDifference_DM_scatter_plot_",
                        filename, ".png"))
    ggsave(p, file=savefile, width=7, height=7)
    
}

check_example <- function(){
    # using this to interactively explore specific runs
    num_sites <- 50 
    num_samples <- 500
    percent_meth_difference <- 0.1
    percent_sites_dm <- 1

    N <- 1000
    dmr_length <- num_sites / 2

    # name of output files are formatted like this
    datadir <- paste0("../../output/01_dmp_results/")
    (filename <- paste0("DM_results",
                        "_numSites", num_sites,
                        "_numSamples", num_samples,
                        "_pct_meth_diff", percent_meth_difference,
                        "_dmr_length", dmr_length,
                        "_pct_sites_dm", percent_sites_dm,
                        "_N", N))
    (datafile <- paste0(datadir, filename, ".rds"))

    results <- readRDS(datafile)
    names(results)

    avgs_res <- results$avgs_res
    names(avgs_res)

    avgs_res$full_meth_dims
    avgs_dm <- avgs_res$results
    head(avgs_dm)
    dim(avgs_dm)

    # little more than 50% identified as DM
    table(avgs_dm$bf_pval < 0.05)
     
    # check the simulations
    datadir <- paste0("../../output/00_simulated_data/")
    (filename <- paste0("simulated_data",
                        "_numSites", num_sites,
                        "_numSamples", num_samples,
                        "_pct_meth_diff", percent_meth_difference,
                        "_dmr_length", dmr_length,
                        "_pct_sites_dm", percent_sites_dm,
                        "_N", N))
    (datafile <- paste0(datadir, filename, ".rds"))

    sims <- readRDS(datafile)
    names(sims)

    sims_dat <- sims$raw_data
    names(sims_dat)

    # check on region
    i <- 1
    simi <- sims_dat[[i]]

    simi_meth <- simi[[1]]
    head(simi_meth)[,1:10]

    # only 1 dmr, drop percent site dm == 1
    simi_dmr <- simi[[2]]
    head(simi_dmr)

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

    # this function looks at percent meth difference and percent sites dm
    # fix these for this analysis
    num_sites <- 20
    num_samples <- 5000
    check_meth_difference(num_sites, num_samples)

    1
}

main()
