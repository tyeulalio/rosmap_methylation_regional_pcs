library(ggplot2)
library(tidyverse)

# example script that processes the simulation results from simulation_test.R script
# this script will look at the results from 1 parameter change in the simulations

savedir <- paste0("../../output/02_processed_results/")
dir.create(savedir)

# use the colors defined in the color panels R script
source("../color_panels.R")

# load the methylation data 
load_data <- function(row){
    percent_meth_difference <- row[['percent_meth_difference']] %>% as.numeric()
    percent_sites_dm <- row[['percent_sites_dm']] %>% as.numeric()
    num_sites <- row[['num_sites']] %>% as.numeric()
    num_samples <- row[['num_samples']] %>% as.numeric()

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

    # add this check for missing files
    # some runs may have failed - will look into these later
    if (!file.exists(datafile)){
        print(filename)
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


    # check the simulations
    testing = FALSE
    if (testing){
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

        # check on region
        i <- 1
        simi <- sims_dat[[i]]

        # methylation values across sites
        simi_meth <- simi[[1]]
        head(simi_meth)[,1:10]

        # number of DM site
        table(simi_meth$dmsites)

        # get counts of DM sites
        dm_sites <- sum(simi_meth$dmsites == 0)
        not_dm_sites <- sum(simi_meth$dm_sites != 0)

        (summarized_res$dm_sites = dm_sites)
        (summarized_res$not_dm_sites = not_dm_sites)
    }

    summarized_res
}

check_parameter <- function(percent_meth_difference_range, percent_sites_dm_range,
                            num_sites_range, num_samples_range, parm, plot){
    runs <- expand.grid(percent_meth_difference=percent_meth_difference_range,
                        percent_sites_dm=percent_sites_dm_range,
                        num_sites=num_sites_range,
                        num_samples=num_samples_range
    ) %>%
        mutate(rownum=row_number()) 
    head(runs)

    # fix these for this analysis
    #num_sites <- 50 
    #num_samples <- 500
    N <- 1000
    dmr_length <- num_sites / 2

    row <- runs[1,] # for testing

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
               num_sites,
               num_samples,
               bh_sig_prop_avgs:bf_sig_prop_rpcs
        ) %>%
        gather('prop_type', 'prop', bh_sig_prop_avgs:bf_sig_prop_rpcs) %>%
        mutate(prop_type = str_remove(prop_type, '_sig_prop')) %>%
        separate(prop_type, c('adjustment', 'summary_type'), sep='_') %>%
        filter(adjustment == 'bf')
    head(long_results)

    if (plot){
        return(long_results)
    }

    parmdir <- paste0(savedir, "parameters_summary_stats/")
    dir.create(parmdir, showWarnings=FALSE)

    save_detailed_results <- function(parm){
        parm_results <- long_results
        parm_results['parm'] <- parm_results[parm]

        # create df to compare averages and rpcs
        detailed_results <- parm_results %>%
            spread(summary_type, prop) %>%
            mutate(diff = rpcs - avgs) %>%
            filter(percent_meth_difference >= 0.01) %>%
            group_by(parm) %>% # change this to parameter we're focusing on
            summarize(min=min(diff),
                      med=median(diff),
                      mean=mean(diff),
                      max=max(diff),
                      avg_min=min(avgs),
                      avg_median=median(avgs),
                      avg_mean=mean(avgs),
                      avg_max=max(avgs),
                      rpc_min=min(rpcs),
                      rpc_median=median(rpcs),
                      rpc_mean=mean(rpcs),
                      rpc_max=max(rpcs)
            )
        head(detailed_results)
        (detailed_results)

        detailed_results$parm_name <- parm

        num_sites <- paste(num_sites_range, collapse='-')
        num_samples  <- paste(num_samples_range, collapse='-')

        savefile <- paste0(parmdir, parm,
                           "_numSites", num_sites,
                           "_numSamples", num_samples,
                           "_simulation_dm_summary_stats.rds")
        saveRDS(detailed_results, savefile)

        detailed_results
    }
    head(long_results)
    res <- save_detailed_results(parm)

    res
}

plot_dm_results <- function(num_sites, num_samples){
    # will vary percent_meth_difference and percent_sites_dm
    #percent_meth_difference_range <- seq(0.1, 0.9, 0.1)
    percent_meth_difference_range <- seq(0.00001, 0.1, 0.01)
    percent_sites_dm_range <- seq(0.25,0.75,0.25)
    #num_sites=50
    #num_samples=50
    num_sites_range=c(num_sites)
    num_samples_range=c(num_samples)
    parm <- 'percent_meth_difference'
    long_results <- check_parameter(percent_meth_difference_range, percent_sites_dm_range, num_sites_range, num_samples_range, parm, plot=TRUE)
    head(long_results)
    dmr_length <- num_sites / 2

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
        #geom_smooth() +
        geom_point() +
        #geom_line() +
        facet_wrap(. ~ percent_sites_dm) +
        #scale_x_continuous(breaks=seq(0, 1, 0.1)) +
        scale_x_continuous(breaks=seq(0, 0.1, 0.02)) +
        ylim(0,1) +
        scale_color_manual(values=summary_type_colors, name="Summary type") + # comes from the color panel file we sourced
        theme_bw() +
        theme(
              text=element_text(size=20),
              axis.text.x=element_text(angle=50, hjust=1, vjust=1),
              # remove minor grid lines
              panel.grid.minor=element_blank(),
              #panel.grid.major=element_blank(), # will do this later for paper
              strip.background=element_blank(), # remove facet strip background
              panel.border=element_rect(color='black',
                                        fill=NA,
                                        linewidth=1
                                        ),
              legend.text=element_text(size=12),
              legend.title=element_text(size=12),
              legend.box.background=element_rect(color='black'),
              legend.position=c(0.87, 0.15)
              ) +
        ylab("Proportion DM detected") +
        xlab("Percent methylation difference")
        #ggtitle(str_replace_all(filename, "_", " "))
    (savefile <- paste0(savedir,
                        #"percentSitesDm_percentMethDifference_DM_scatter_plot_",
                        "smallRange_percentSitesDm_percentMethDifference_DM_scatter_plot_",
                        #"verySmallRange_percentSitesDm_percentMethDifference_DM_scatter_plot_",
                        filename, ".png"))
    ggsave(p, file=savefile, width=7, height=5)
    

}

save_parameter_summary_stats <- function(num_sites, num_samples){
    # will vary percent_meth_difference and percent_sites_dm
    #percent_meth_difference_range <- seq(0.1, 0.9, 0.1)
    percent_meth_difference_range <- seq(0.00001, 0.1, 0.01)
    percent_sites_dm_range <- seq(0.25,0.75,0.25)
    #num_sites=50
    #num_samples=50
    results <- list(test=0)

    num_sites_range=c(num_sites)
    num_samples_range=c(num_samples)
    parm <- 'percent_meth_difference'
    results[[parm]] <- check_parameter(percent_meth_difference_range, percent_sites_dm_range, num_sites_range, num_samples_range, parm, plot=FALSE)

    parm <- 'percent_sites_dm'
    results[[parm]] <- check_parameter(percent_meth_difference_range, percent_sites_dm_range, num_sites_range, num_samples_range, parm, plot=FALSE)

    num_sites_range=c(20,50)
    parm <- 'num_sites'
    results[[parm]] <- check_parameter(percent_meth_difference_range, percent_sites_dm_range, num_sites_range, num_samples_range, parm, plot=FALSE)

    num_sites_range=c(num_sites)
    num_samples_range=c(50,500)
    parm <- 'num_samples'
    results[[parm]] <- check_parameter(percent_meth_difference_range, percent_sites_dm_range, num_sites_range, num_samples_range, parm, plot=FALSE)


    parmres_df <- do.call(rbind, results)[-1,] %>%
        mutate(num_sites=num_sites, num_samples=num_samples) %>%
        select(parameter_name=parm_name, parameter_value=parm, 
               num_sites, num_samples,
               diff_min=min, diff_med=med, diff_mean=mean, diff_max=max,
               everything()) 
    head(parmres_df)
    dim(parmres_df)


    savefile <- paste0(parmdir, "numSites", num_sites, "_numSamples", num_samples,
                       "_allParameters_summary_stats.csv"
    )
    savefile
    write_csv(parmres_df, savefile)

    1
}

check_example <- function(){
    # using this to interactively explore specific runs
    num_sites <- 50 
    num_samples <- 500
    percent_sites_dm <- 0.75

    (percent_meth_difference_range <- seq(0.00001, 0.006, 0.0006))
    (percent_meth_difference <- percent_meth_difference_range[[5]])

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
    rpcs_dm <- results$rpcs_res$results
    head(avgs_dm)
    head(rpcs_dm)

    # little more than 50% identified as DM
    table(avgs_dm$bf_pval < 0.05)
    table(results$rpcs_res$results$bf_pval < 0.05)
    quantile(results$avgs_res$results$P.Value, seq(0, 1, 0.1))
    quantile(results$rpcs_res$results$P.Value, seq(0, 1, 0.1))
     
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

    # check on region
    i <- 1
    simi <- sims_dat[[i]]

    # methylation values across sites
    simi_meth <- simi[[1]]
    head(simi_meth)[,1:10]

    # number of DM site
    table(simi_meth$dmsites)

    # only 1 dmr, drop percent site dm == 1
    simi_dmr <- simi[[2]]
    head(simi_dmr)

}

combine_summary_stats <- function(){
    runs <- expand.grid(num_sites=c(20,50),
                        num_samples=c(50,500)
    )
    load_stats <- function(run){
        num_sites <- run[['num_sites']] %>% as.numeric()
        num_samples <- run[['num_samples']] %>% as.numeric()

        parmdir <- paste0(savedir, "parameters_summary_stats/")
        savefile <- paste0(parmdir, "numSites", num_sites, "_numSamples", num_samples,
                           "_allParameters_summary_stats.csv"
        )
        parmres_df <- read_csv(savefile)
    }
    res <- apply(runs, 1, load_stats)
    res_df <- do.call(rbind, res)
    head(res_df)
    dim(res_df)

    savefile <- paste0(parmdir, "simulation_summary_stats.csv")
    savefile
    write_csv(res_df, savefile)
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
    N <- 1000
    num_sites <- 50
    num_samples <- 500
    parmdir <- paste0(savedir, "parameters_summary_stats/")
    save_parameter_summary_stats(num_sites, num_samples)
    plot_dm_results(num_sites, num_samples)


    1
}

main()
