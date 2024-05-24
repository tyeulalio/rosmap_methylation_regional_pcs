library(ggplot2)
library(PRROC)
library(GenomicRanges)
library(tidyverse)

# example script that processes the simulation results from simulation_test.R script
# this script will look at the results from 1 parameter change in the simulations

savedir <- paste0("../../output/04_simulation_performance_metrics/")
dir.create(savedir)

# use the colors defined in the color panels R script
source("../color_panels.R")

load_dmr_data <- function(runname){
    # load dmr-level data
    datadir <- paste0("../../output/03_simulation_dmrs/")
    list.files(datadir) %>% head()

    (datafile <- paste0(datadir,
                       "dmrcate_results",
                       runname,
                       ".rds"
    ))
    # add this check for missing files
    # some runs may have failed - will look into these later
    if (!file.exists(datafile)){
        print(datafile)
        print("datafile does not exist")
        return(NA)
    }
    dm_res <- readRDS(datafile) 
    head(dm_res)

    dm_res
}

load_cpg_data <- function(runname){
    # load cpg-level data
    datadir <- paste0("../../output/03_simulation_dmrs/")
    #list.files(datadir)

    (datafile <- paste0(datadir,
                       "limma_dmp_results",
                       runname,
                       ".rds"
    ))
    # add this check for missing files
    # some runs may have failed - will look into these later
    if (!file.exists(datafile)){
        print(filename)
        print("datafile does not exist")
        return(NA)
    }
    dm_res <- readRDS(datafile)
    head(dm_res)

    dm_res
}

# load the methylation data 
load_dm_results <- function(runname){
    # name of output files are formatted like this
    datadir <- paste0("../../output/01_dmp_results/")
    (filename <- paste0("DM_results", runname))
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

    # get the DM results
    rpcs_res <- results[['rpcs_res']]$results %>%
        mutate(stat = bh_pval)
    avgs_res <- results[['avgs_res']]$results %>%
        mutate(stat = bh_pval)
    head(rpcs_res)
    head(avgs_res)

    # load the cpg-level and dmr-level data
    cpg_res <- load_cpg_data(runname) %>%
        mutate(stat = bh_pval)
    head(cpg_res)

    # load dmr results
    dmr_res <- load_dmr_data(runname) 
    if (!all(is.na(dmr_res))){
        dmr_res <- dmr_res %>%
            mutate(stat = HMFDR)
    }
    head(dmr_res)

    dm_res <- list(rpcs=rpcs_res, avgs=avgs_res, cpgs=cpg_res, dmr=dmr_res)
    return(dm_res)

    # check the simulations
    testing = FALSE
    if (testing){
        datadir <- paste0("../../output/00_simulated_data/")
        (filename <- paste0("simulated_data", runname))
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

    dm_res
}

load_negative_dm_results <- function(runname){
    # name of output files are formatted like this
    datadir <- paste0("../../output/01_dmp_results/")
    (filename <- paste0("randomizedPheno_DM_results", runname))
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

    # get the DM results
    rpcs_res <- results[['rpcs_res']]$results
    avgs_res <- results[['avgs_res']]$results
    head(rpcs_res)
    head(avgs_res)

    dm_res <- list(rpcs=rpcs_res, avgs=avgs_res)

    dm_res
}

format_data <- function(positive_dm_res, negative_dm_res){
    # format the results for pr curves
    names(positive_dm_res)

    # add label to the data
    pos_res <- lapply(positive_dm_res, function(x) x %>% mutate(label = 1))
    neg_res <- lapply(negative_dm_res, function(x) x %>% mutate(label = 0))

    # format data for each summary type
    st <- 'rpcs'
    get_scores <- function(pos_res, neg_res, st) {
        # combine positive and negative results
        combined_res <- rbind(pos_res[[st]], neg_res[[st]]) %>%
            # create score for pr curve
            mutate(score = -log10(bh_pval)) 
        combined_res
    }
    rpcs_res <- get_scores(pos_res, neg_res, 'rpcs')
    avgs_res <- get_scores(pos_res, neg_res, 'avgs')
    cpgs_res <- get_scores(pos_res, neg_res, 'cpgs')

    head(avgs_res)
    head(cpgs_res)

    # make sure thre's only one score per region for pcs and cpgs
    cpgs <- cpgs_res %>%
        rownames_to_column('region_cpg') %>%
        separate(region_cpg, c('region', 'cpg'), sep='_') %>%
        select(region, score) %>%
        group_by(region) %>%
        summarize(score = max(score)) %>%
        mutate(label = cpgs_res[1,'label'])
    head(cpgs)
    cpgs_res <- cpgs

    head(rpcs_res)
    rpcs <- rpcs_res %>%
        rownames_to_column('region_pc') %>%
        separate(region_pc, c('region', 'pc'), sep='-') %>%
        select(region, score) %>%
        group_by(region) %>%
        summarize(score = max(score)) %>%
        mutate(label = rpcs_res[1,'label'])
    head(rpcs)
    rpcs_res <- rpcs

    pr_data <- list(avgs=avgs_res, rpcs=rpcs_res, cpgs=cpgs_res)
    pr_data
}

# get precision recall curve
#dm_res <- roc_data[[st]]
compute_roc <- function(dm_res, summary_type){
    # create precision recall curve
    roc <- roc.curve(scores.class0 = dm_res$score[dm_res$label == 1],
                   scores.class1 = dm_res$score[dm_res$label == 0],
                   curve = TRUE
    )
    names(roc)
    head(roc$curve)

    # extract necessary data
    roc_curve <- data.frame(false_pos_rate = roc$curve[,1],
                          sensitivity = roc$curve[,2]
    ) %>%
        mutate(summary_type = summary_type) 
    head(roc_curve)

    roc[['roc_curve']] <- roc_curve
    roc
}

roc_analysis <- function(positive_dm_res, negative_dm_res, runname){
    # plot precision recall curves

    # format the data
    roc_data <- format_data(positive_dm_res, negative_dm_res)


    st <- 'avgs'
    avgs <- compute_roc(roc_data[[st]], st) 
    st <- 'rpcs'
    rpcs <- compute_roc(roc_data[[st]], st) 
    st <- 'cpgs'
    cpgs <- compute_roc(roc_data[[st]], st) 

    names(avgs)

    combined_roc <- list(avgs=avgs, rpcs=rpcs, cpgs=cpgs)

    (savefile <- paste0(savedir,
                       "ROC_results",
                       runname,
                       ".rds"
    ))
    saveRDS(combined_roc, savefile)

    # plot the pr curve 
    #prcurves <- rbind(avgs_pr$pr_curve, rpcs_pr$pr_curve)
    #head(prcurves)
    #plot_pr_curve(prcurves)

    1
}



# get precision recall curve
#dm_res <- pr_data[[st]]
compute_pr <- function(dm_res, summary_type){
    # create precision recall curve
    pr <- pr.curve(scores.class0 = dm_res$score[dm_res$label == 1],
                   scores.class1 = dm_res$score[dm_res$label == 0],
                   curve = TRUE
    )
    names(pr)
    head(pr$curve)

    # extract necessary data
    pr_curve <- data.frame(recall = pr$curve[,1],
                          precision = pr$curve[,2]
    ) %>%
        mutate(summary_type = summary_type) 
    head(pr_curve)

    pr[['pr_curve']] <- pr_curve
    pr
}

plot_pr_curve <- function(prcurves){
    # plot the pr curves 
    p <- ggplot(prcurves) +
        geom_line(aes(x=recall,
                      y=precision,
                      color=summary_type)) +
        labs(title=paste0("Precision-Recall Curve"),
             x = "Recall",
             y = "Precision"
        ) +
    theme_bw()

    (savefile <- paste0(savedir, 
                       "pr_curve",
                       runname, 
                       ".png"
                    ))
    ggsave(p, file=savefile)
}


pr_analysis <- function(positive_dm_res, negative_dm_res, runname){
    # plot precision recall curves

    # format the data
    pr_data <- format_data(positive_dm_res, negative_dm_res)

    st <- 'avgs'
    avgs_pr <- compute_pr(pr_data[[st]], st) 
    st <- 'rpcs'
    rpcs_pr <- compute_pr(pr_data[[st]], st) 
    st <- 'cpgs'
    cpgs_pr <- compute_pr(pr_data[[st]], st) 

    names(avgs_pr)

    combined_pr <- list(avgs=avgs_pr, rpcs=rpcs_pr, cpgs=cpgs_pr)

    (savefile <- paste0(savedir,
                       "PR_results",
                       runname,
                       ".rds"
    ))
    saveRDS(combined_pr, savefile)

    # plot the pr curve 
    #prcurves <- rbind(avgs_pr$pr_curve, rpcs_pr$pr_curve)
    #head(prcurves)
    #plot_pr_curve(prcurves)

    1
}

plot_aurocs <- function(num_sites, num_samples){
    # plot the auprcs together

    percent_meth_difference_range <- c(seq(0.00001, 0.1, 0.01))
    percent_sites_dm_range <- seq(0.25,0.75,0.25)

    #num_sites <- 50
    #num_samples <- 500

    num_sites_range=c(num_sites)
    #num_sites_range=c(20, 50)
    num_samples_range=c(num_samples)
    #num_samples_range=c(50,500)
    #num_samples_range=c(50,500,5000)

    runs <- expand.grid(percent_meth_difference=percent_meth_difference_range,
                        percent_sites_dm=percent_sites_dm_range,
                        num_sites=num_sites_range,
                        num_samples=num_samples_range
    ) %>%
        mutate(rownum=row_number()) 
    head(runs)
    dim(runs)

    # need to load the auprc for each run
    run <- run[1,]
    load_auroc <- function(run){
        runname <- get_runname(run)

        savefile <- paste0(savedir,
                           "ROC_results",
                           runname,
                           ".rds"
        )
        roc_res <- readRDS(savefile)
        names(roc_res$avgs)

        # select auc measure
        auc_measure <- 'auc'

        # get aucs
        aucs <- data.frame(summary_type=c('avgs', 'rpcs', 'cpgs'),
                           auc=c(roc_res[['avgs']][[auc_measure]], 
                                 roc_res[['rpcs']][[auc_measure]],
                                 roc_res[['cpgs']][[auc_measure]]
                           )
        ) 
        head(aucs)

        run_parameters <- get_run_parameters(run)
        names(run_parameters)

        aucs <- aucs %>%
            mutate(percent_meth_difference = run_parameters[['percent_meth_difference']],
                   percent_sites_dm = run_parameters[['percent_sites_dm']],
                   num_sites = run_parameters[['num_sites']],
                   num_samples = run_parameters[['num_sampes']],
                   dmr_length = run_parameters[['dmr_length']]
            )
        head(aucs)

        aucs
    }

    aucs <- apply(runs, 1, load_auroc)
    aucs_df <- do.call(rbind, aucs)
    head(aucs_df)
    
    # plot aucs
    p <- ggplot(aucs_df) +
        geom_point(aes(x=percent_meth_difference,
                       y=auc,
                       color=summary_type
                       )) +
        facet_wrap(. ~ percent_sites_dm) +
        labs(title = paste("AUROC -", 
                           num_sites, "sites,",
                           num_samples, "samples"
                           ),
             y = "AUROC", 
             x = "Percent methylation difference"
             ) +
        theme_bw()

    (savefile <- paste0(savedir, 
                       "AUROC_scatter_plot",
                       "_numSites", num_sites,
                       "_numSamples", num_samples,
                       ".png"
    ))
    ggsave(p, file=savefile, height=6, width=7)
}

plot_auprcs <- function(num_sites, num_samples){
    # plot the auprcs together

    percent_meth_difference_range <- c(seq(0.00001, 0.1, 0.01))
    percent_sites_dm_range <- seq(0.25,0.75,0.25)

    #num_sites <- 50
    #num_samples <- 500

    num_sites_range=c(num_sites)
    #num_sites_range=c(20, 50)
    num_samples_range=c(num_samples)
    #num_samples_range=c(50,500)
    #num_samples_range=c(50,500,5000)

    runs <- expand.grid(percent_meth_difference=percent_meth_difference_range,
                        percent_sites_dm=percent_sites_dm_range,
                        num_sites=num_sites_range,
                        num_samples=num_samples_range
    ) %>%
        mutate(rownum=row_number()) 
    head(runs)
    dim(runs)

    # need to load the auprc for each run
    run <- runs[1,]
    load_auprc <- function(run){
        runname <- get_runname(run)

        savefile <- paste0(savedir,
                           "PR_results",
                           runname,
                           ".rds"
        )
        pr_res <- readRDS(savefile)
        names(pr_res$avgs)

        # select auc measure
        auc_measure <- 'auc.integral'

        # get aucs
        aucs <- data.frame(summary_type=c('avgs', 'rpcs', 'cpgs'),
                           auc=c(pr_res[['avgs']][[auc_measure]], 
                                 pr_res[['rpcs']][[auc_measure]],
                                 pr_res[['cpgs']][[auc_measure]]
                           )
        ) 
        head(aucs)

        run_parameters <- get_run_parameters(run)
        names(run_parameters)

        aucs <- aucs %>%
            mutate(percent_meth_difference = run_parameters[['percent_meth_difference']],
                   percent_sites_dm = run_parameters[['percent_sites_dm']],
                   num_sites = run_parameters[['num_sites']],
                   num_samples = run_parameters[['num_sampes']],
                   dmr_length = run_parameters[['dmr_length']]
            )
        head(aucs)

        aucs
    }

    aucs <- apply(runs, 1, load_auprc)
    aucs_df <- do.call(rbind, aucs)
    head(aucs_df)
    
    # plot aucs
    p <- ggplot(aucs_df) +
        geom_point(aes(x=percent_meth_difference,
                       y=auc,
                       color=summary_type
                       )) +
        facet_wrap(. ~ percent_sites_dm) +
        labs(title = paste("AUPRC -", 
                           num_sites, "sites,",
                           num_samples, "samples"
                           ),
             y = "AUPRC", 
             x = "Percent methylation difference"
             ) +
        theme_bw()

    (savefile <- paste0(savedir, 
                       "AUPRC_scatter_plot",
                       "_numSites", num_sites,
                       "_numSamples", num_samples,
                       ".png"
    ))
    ggsave(p, file=savefile, height=6, width=7)
}

get_run_parameters <- function(run){
    percent_meth_difference <- run[['percent_meth_difference']] %>% as.numeric()
    percent_sites_dm <- run[['percent_sites_dm']] %>% as.numeric()
    num_sites <- run[['num_sites']] %>% as.numeric()
    num_samples <- run[['num_samples']] %>% as.numeric()

    dmr_length <- num_sites / 2
    N <- 1000

    list(percent_meth_difference=percent_meth_difference,
         percent_sites_dm=percent_sites_dm,
         num_sites=num_sites,
         num_samples=num_samples,
         dmr_length=dmr_length,
         N=N
    )
}

get_runname <- function(run){
    percent_meth_difference <- run[['percent_meth_difference']] %>% as.numeric()
    percent_sites_dm <- run[['percent_sites_dm']] %>% as.numeric()
    num_sites <- run[['num_sites']] %>% as.numeric()
    num_samples <- run[['num_samples']] %>% as.numeric()


    dmr_length <- num_sites / 2
    N <- 1000

    # name of output files are formatted like this
    (runname <- paste0(
                        "_numSites", num_sites,
                        "_numSamples", num_samples,
                        "_pct_meth_diff", percent_meth_difference,
                        "_dmr_length", dmr_length,
                        "_pct_sites_dm", percent_sites_dm,
                        "_N", N))
    runname
}

# compute metrics
compute_metrics <- function(dm_res, params){
    # compute TP, FP, TN, FN results
    names(params)
    head(dm_res$rpcs)

    sig_thresh <- 0.05

    # get num sig and total for rpcs
    rpcs <- dm_res$rpcs %>%
        rownames_to_column('region_pc') %>%
        separate(region_pc, c('region', 'pc'), sep='-') %>%
        select(region, stat) %>%
        group_by(region) %>%
        summarize(stat = min(stat)) %>%
        unique()
    head(rpcs)
    rpcs_counts <- data.frame(summary_type='rpcs', 
                              num_sig = sum(rpcs$stat < sig_thresh),
                              total_regions = length(unique(rpcs$region))
    )
    rpcs_counts

    # get num sig and total for avgs
    avgs <- dm_res$avgs 
    head(avgs)
    avgs_counts <- data.frame(summary_type='avgs',
                              num_sig = sum(avgs$stat < sig_thresh),
                              total_regions = nrow(avgs)
    )
    avgs_counts

    # get num sig and total for cpgs
    cpgs <- dm_res$cpgs %>%
        rownames_to_column('region_cpg') %>%
        separate(region_cpg, c('region', 'cpg'), sep="_") %>%
        select(region, stat) %>%
        group_by(region) %>%
        summarize(stat = min(stat)) %>%
        unique()
    head(cpgs)
    cpgs_counts <- data.frame(summary_type='cpgs',
                              num_sig = sum(cpgs$stat < sig_thresh),
                              total_regions = nrow(cpgs)
    )
    cpgs_counts

    # combine the counts
    all_counts <- rpcs_counts %>%
        full_join(avgs_counts) %>%
        full_join(cpgs_counts) 

    # get num sig and total for dmrs
    if (all(is.na(dm_res$dmr))){
        NA
    } else{
        dmrs <- dm_res$dmr %>%
            separate(coord, c('region', 'dmr'), sep=':') %>%
            select(region, stat) %>%
            group_by(region) %>%
            summarize(stat = min(stat)) %>%
            unique()
        head(dmrs)
        dmrs_counts <- data.frame(summary_type = 'dmrs',
                                  num_sig = sum(dmrs$stat < sig_thresh),
                                  total_regions = nrow(dmrs)
        )
        dmrs_counts

        # combine with other results
        all_counts <- all_counts %>%
            full_join(dmrs_counts)
    }

    all_counts
    
    if (params$percent_meth_difference == 0){
        metrics <- all_counts %>%
            mutate(tp = 0,
                   fp = num_sig,
                   tn = total_regions - num_sig,
                   fn = 0
            )
    } else{
        metrics <- all_counts %>%
            mutate(tp = num_sig,
                   fp = 0,
                   tn = 0,
                   fn = total_regions - num_sig
            )
    }
    metrics
}

#dm_res <- positive_dm_res
compare_dm_results <- function(dm_res, run){
    # compare the DM results
    params <- get_run_parameters(run)
    params

    metrics <- compute_metrics(dm_res, params)
    metrics

    names(params)

    # add parameters into the data frame
    metrics_params <- metrics %>%
        mutate(percent_meth_difference = params$percent_meth_difference,
               percent_sites_dm = params$percent_sites_dm,
               num_sites = params$num_sites,
               num_samples = params$num_samples,
               dmr_length = params$dmr_length
        )
    metrics_params

    runname <- get_runname(run)
    (savefile <- paste0(savedir,
                        "dm_summary_metrics",
                        runname,
                        ".rds"
                        ))
    saveRDS(metrics_params, savefile)
    1
}

plot_metrics <- function(num_sites, num_samples){
    # plot summary metrics
    #percent_meth_difference_range <- c(0, seq(0.00001, 0.1, 0.01))
    percent_meth_difference_range <- c(0, seq(0.01001, 0.1, 0.01))
    percent_sites_dm_range <- seq(0.25,0.75,0.25)

    #num_sites <- 50
    #num_samples <- 500

    num_sites_range=c(num_sites)
    #num_sites_range=c(20, 50)
    num_samples_range=c(num_samples)
    #num_samples_range=c(50,500)
    #num_samples_range=c(50,500,5000)

    runs <- expand.grid(percent_meth_difference=percent_meth_difference_range,
                        percent_sites_dm=percent_sites_dm_range,
                        num_sites=num_sites_range,
                        num_samples=num_samples_range
    ) %>%
        mutate(rownum=row_number()) 
    head(runs)
    dim(runs)

    # need to load the auprc for each run
    run <- run[1,]
    process_run <- function(run){
        runname <- get_runname(run)
        (datafile <- paste0(savedir,
                            "dm_summary_metrics",
                            runname,
                            ".rds"
                            ))
        metrics <- readRDS(datafile)
        metrics
    }

    metrics <- apply(runs, 1, process_run)
    metrics_df <- do.call(rbind, metrics)
    head(metrics_df)

    metrics_df <- metrics_df %>%
        mutate(prop_df = num_sig / total_regions)
    head(metrics_df)

    metrics_df %>%
        select(summary_type, percent_meth_difference:num_samples, prop_df) %>%
        spread(summary_type, prop_df)

    # plot the metrics
    p <- ggplot(metrics_df) +
        geom_point(aes(x=percent_meth_difference,
                       y=prop_df,
                       color=summary_type)) +
        facet_wrap(. ~ percent_sites_dm) +
        labs(title=paste(num_sites, "sites,", num_samples, "samples"),
             x="Percent methylation difference"
             #y="Proportion DM detected"
        ) +
        theme_bw()

    (savefile <- paste0(savedir,
                        "proportion_DM_scatter",
                       "_numSites", num_sites,
                       "_numSamples", num_samples,
                       ".png"
    ))
    ggsave(p, file=savefile)
    
    1
}

get_true_dmr_cpgs <- function(run){
    # get positions of the true dmrs from simulations
    runname <- get_runname(run)

    datadir <- paste0("../../output/00_simulated_data/")
    (datafile <- paste0(datadir,
                       "simulated_data",
                       runname,
                       ".rds"
    ))
    simulated_data <- readRDS(datafile)
    names(simulated_data)

    # get raw simulated data
    raw_data <- simulated_data[['raw_data']]
    head(raw_data[[1]][[1]])[1:10]

    # need to grab the second dataframe from each list item
    i <- 1000
    get_dmr_sites <- function(i){
        #print(i)
        item <- raw_data[[i]]
        if (all(is.na(item))) return(NA)

        # grab cpg labels
        head(item[[1]])[,1:10]
        cpg_labels <- item[[1]] %>%
            select(pos=chrpos, dm=dmsites) %>%
            mutate(region=paste0('gene', i))
        head(cpg_labels)

        #dmr_sites <- raw_data[[i]][[2]] %>%
            #as.data.frame()

        #if (all(is.na(dmr_sites))) return(NA)

        ## case where there is only 1 df2 row
        #if (ncol(dmr_sites) == 1){
            #dmr_sites <- dmr_sites %>%
                #t() %>%
                #as.data.frame()
            #names(dmr_sites) <- c('dmr_row_start', 'dmr_chrpos', 'dmr_number_sites')
        #}
        #dmr_sites

        #dmr_sites <- dmr_sites %>%
            #mutate(region = i)
        #dmr_sites

        cpg_labels
    }
    dmr_sites <- lapply(1:length(raw_data), get_dmr_sites)
    #full_dmr_sites <- dmr_sites[!all(is.na(dmr_sites))]
    dmr_sites_df <- do.call(rbind, dmr_sites)
    head(dmr_sites_df)

    dmr_sites_df
}


get_true_dmr_cpgs <- function(run){
    # get positions of the true dmrs from simulations
    runname <- get_runname(run)

    datadir <- paste0("../../output/00_simulated_data/")
    (datafile <- paste0(datadir,
                       "simulated_data",
                       runname,
                       ".rds"
    ))
    simulated_data <- readRDS(datafile)
    names(simulated_data)

    # get raw simulated data
    raw_data <- simulated_data[['raw_data']]
    head(raw_data[[1]][[1]])[1:10]

    # need to grab the second dataframe from each list item
    i <- 1000
    get_dmr_sites <- function(i){
        #print(i)
        item <- raw_data[[i]]
        if (all(is.na(item))) return(NA)

        # grab cpg labels
        head(item[[1]])[,1:10]
        cpg_labels <- item[[1]] %>%
            select(pos=chrpos, dm=dmsites) %>%
            mutate(region=paste0('gene', i))
        head(cpg_labels)

        #dmr_sites <- raw_data[[i]][[2]] %>%
            #as.data.frame()

        #if (all(is.na(dmr_sites))) return(NA)

        ## case where there is only 1 df2 row
        #if (ncol(dmr_sites) == 1){
            #dmr_sites <- dmr_sites %>%
                #t() %>%
                #as.data.frame()
            #names(dmr_sites) <- c('dmr_row_start', 'dmr_chrpos', 'dmr_number_sites')
        #}
        #dmr_sites

        #dmr_sites <- dmr_sites %>%
            #mutate(region = i)
        #dmr_sites

        cpg_labels
    }
    dmr_sites <- lapply(1:length(raw_data), get_dmr_sites)
    #full_dmr_sites <- dmr_sites[!all(is.na(dmr_sites))]
    dmr_sites_df <- do.call(rbind, dmr_sites)
    head(dmr_sites_df)

    dmr_sites_df
}

get_cpg_labels <- function(cpgs, dmr_label){
    # add cpg-level labels
    head(cpgs)
    cpg_labels <- cpgs %>%
        rownames_to_column('region_cpg') %>%
        separate(region_cpg, c('region', 'cpg'), sep='_', remove=FALSE) %>%
        left_join(dmr_label) %>%
        select(-region, -cpg) %>%
        column_to_rownames('region_cpg')
    head(cpg_labels)
    cpg_labels
}

get_dmr_labels <- function(dmrs, dmr_label){
    # add dmr labels
    head(dmrs)
    head(dmr_label)

    # need to find overlap between dmrs and true dmrs
    formatted_dmrs <- dmrs %>%
        separate(coord, c('chr', 'dmr'), sep=':', remove=FALSE) %>%
        separate(dmr, c('start', 'end'), sep='-')
    head(formatted_dmrs)

    head(dmr_label)
    formatted_true <- dmr_label %>%
        rename(chr = region) %>%
        mutate(start = cpg, end = cpg)
    head(formatted_true)

    # create genomic ranges objects
    dmrs_gr <- makeGRangesFromDataFrame(formatted_dmrs, keep.extra.columns = TRUE)
    true_gr <- makeGRangesFromDataFrame(formatted_true, keep.extra.columns = TRUE)

    # find overlaps
    overlaps <- findOverlaps(subject=dmrs_gr, query=true_gr) %>%
        as.data.frame()
    head(overlaps)

    dmr_cpg_labels <- cbind(as.data.frame(dmrs_gr)[overlaps$subjectHits,],
                         true=as.data.frame(true_gr)[overlaps$queryHits,]
    ) %>%
        select(coord, cpg_pos=true.start, cpg_label=true.cpg_label) 
    head(dmr_cpg_labels)

    dmrs_labelled <- dmrs %>%
        left_join(dmr_cpg_labels)
    head(dmrs_labelled)     

    # check how many cpgs per region are actually DM
    head(dmrs_labelled)
    cpg_counts <- dmrs_labelled %>%
        group_by(coord) %>%
        summarize(true_dm_cpgs=sum(cpg_label),
                  dmr_cpgs=unique(no.cpgs)
        ) %>%
        mutate(prop_true_dm = true_dm_cpgs/dmr_cpgs)
    head(cpg_counts)
    summary(cpg_counts$prop_true_dm)

    # are there any true dm cpgs that are not in DMRs?
    head(formatted_true)
    head(dmr_cpg_labels) # overlaps between dmrs and cpgs
    true_missing <- formatted_true %>%
        mutate(true_cpg_label = cpg_label,
               cpg_pos = cpg
               ) %>%
        select(cpg_pos, true_cpg_label, region=chr) 
    head(true_missing)

    dmr_cpgs <- dmr_cpg_labels %>%
        separate(coord, c('region', 'dmr'), sep=':')
    head(dmr_cpgs)

    true_cpgs_dmrs <- true_missing %>%
        mutate(cpg_pos = as.numeric(cpg_pos)) %>%
        left_join(dmr_cpgs)
    head(true_cpgs_dmrs)

    missing_true_cpgs <- true_cpgs_dmrs %>%
        filter(is.na(dmr),
               true_cpg_label == 1
        )
    head(missing_true_cpgs)
    nrow(missing_true_cpgs)

    dim(dmrs_labelled)
    dim(dmrs)

}

#dm_res <- positive_dm_res
add_label <- function(dm_res, run){
    # add postive/negative label to results
    names(dm_res)

    # rpcs and avgs labels always depends on the run
    params <- get_run_parameters(run)
    names(params)

    # label = 1 if percent meth difference is not zero
    (region_label <- as.numeric(params[['percent_meth_difference']] != 0))
    avgs <- dm_res[['avgs']]  %>%
        mutate(region_label = region_label)
    head(avgs)
    rpcs <- dm_res[['rpcs']]  %>%
        mutate(region_label = region_label)
    head(rpcs)

    # add the region label to the cpgs and dmr regions
    cpgs <- dm_res[['cpgs']] %>%
        mutate(region_label = region_label)
    head(cpgs)

    dmrs <- dm_res[['dmr']]
    if (!is.null(dmrs)){
        dmrs <- dmrs %>%
            mutate(region_label = region_label)
    }
    head(dmrs)

    # add cpg-level labels
    true_dmrs <- get_true_dmr_cpgs(run)
    head(true_dmrs)

    dmr_label <- true_dmrs %>%
        mutate(cpg = as.character(pos),
               region) %>%
        select(region, cpg, cpg_label=dm) %>%
        mutate(cpg_label = abs(cpg_label))
    head(dmr_label)

    cpgs <- get_cpg_labels(cpgs, dmr_label)
    head(cpgs)

    # add labels to dmrs
    dmr_labels <- get_dmr_labels(dmrs, dmr_label)
    head(dmr_labels)
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


    #percent_meth_difference_range <- seq(0.1, 0.9, 0.1)
    percent_meth_difference_range <- c(seq(0.00001, 0.1, 0.01))
    percent_sites_dm_range <- seq(0.25,0.75,0.25)
    num_sites_range=c(20, 50)
    num_samples_range=c(50,500)
    #num_samples_range=c(50,500,5000)

    runs <- expand.grid(percent_meth_difference=percent_meth_difference_range,
                        percent_sites_dm=percent_sites_dm_range,
                        num_sites=num_sites_range,
                        num_samples=num_samples_range
    ) %>%
        mutate(rownum=row_number()) 
    head(runs)
    dim(runs)

    run <- runs[60,]
    process_run <- function(run){
        print(paste("Processing run", run[['rownum']], "out of", nrow(runs)))

        runname <- get_runname(run)

        # load the DM reults
        positive_dm_res <- load_dm_results(runname)
        names(positive_dm_res)
        #head(positive_dm_res[['dmr']])

        # load negative dm results
        negative_run <- run
        negative_run['percent_meth_difference'] = 0
        negative_runname <- get_runname(negative_run)
        negative_dm_res <- load_dm_results(negative_runname)
        #negative_dm_res <- load_negative_dm_results(runname)

        # add labels
        add_label(positive_dm_res, run)

        # compute and save summary metrics
        compare_dm_results(positive_dm_res, run)
        compare_dm_results(negative_dm_res, negative_run)

        # remove the dmrs since we don't have sig values
        # for negative results
        positive_dm_res <- positive_dm_res[1:3]
        negative_dm_res <- negative_dm_res[1:3]

        pr_analysis(positive_dm_res, negative_dm_res, runname)
        roc_analysis(positive_dm_res, negative_dm_res, runname)

        1 
    }
    res <- apply(runs, 1, process_run)

    # plot out AUPRCs
    #num_sites_range=c(20, 50)
    #num_samples_range=c(50,500,5000)
    num_sites <- 20
    num_samples <- 50
    plot_auprcs(num_sites, num_samples)
    plot_aurocs(num_sites, num_samples)

    plot_metrics(num_sites, num_samples)

    1
}

main()
