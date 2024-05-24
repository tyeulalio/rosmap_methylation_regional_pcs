# Run differential methylation analysis on the simulated regions

library(RNOmni)
library(ggplot2)
library(limma)
library(tidyverse)


# create a directory to save the data
savedir <- paste0("../../output/01_dmp_results/")
dir.create(savedir)

datadir <- paste0("../../output/00_simulated_data/")

source("../color_panels.R")


#meth <- rpcs # for testing
run_dmp <- function(meth, pheno, transform=TRUE){
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
    print("applying inverse normal transform")
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

run_leveled_mixed_dmp <- function(){
    # combine DM runs with no-DM runs to test sensitivity/specificity

    # select parameters
    percent_meth_low <- 0

    num_sites <- 50
    num_samples <- 500
    percent_sites_dm <- 0.5
    dmr_length <- num_sites / 2
    N <- 1000

    percent_meth_difference <- 0
    load_simulated_data <- function(percent_meth_difference){
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
        
        names(results)

        # get labels for regions about DM
        simulations <- results$raw_data
        length(simulations)

        i <- 1
        check_region_dm_status <- function(i){
            simi <- simulations[[i]]
            if (all(is.na(simi))) return(NA)

            simi[[2]]

            # get the DM status of sites in this region
            dm_sites <- simi[[1]][,1:2]
            head(dm_sites)

            # count how many sites were DM
            tibble(num_dm_sites=sum(dm_sites$dmsites != 0),
                 num_sites_total=nrow(dm_sites))

        }
        dm_sites_list <- lapply(1:length(simulations), check_region_dm_status)
        dm_sites <- do.call(rbind, dm_sites_list) %>%
            mutate(region = row_number()) 
        head(dm_sites)

        summary(dm_sites$num_dm_sites)


        # grab the averages and rpcs
        avgs <- results$avgs
        rpcs <- results$rpcs

        rownames(avgs) <- paste(rownames(avgs), percent_meth_difference, sep='_')
        rownames(rpcs) <- paste(rownames(rpcs), percent_meth_difference, sep='_')

        rpcs = na.omit(rpcs)
        avgs = na.omit(avgs)

        list(avgs=avgs, rpcs=rpcs, dm_sites=dm_sites)
    }
    
    process_run <- function(){
        dmr_length <- num_sites / 2
        N <- 1000

        # generate a mixed data frame
        no_dm_count  <- N * 0.7
        dm_count <- N * 0.3

        # get simulated data 
        no_dm_data <- load_simulated_data(0)

        # get the no dm summary types
        no_dm_avgs <- no_dm_data$avgs[1:no_dm_count,]
        no_dm_rpcs <- no_dm_data$rpcs[1:no_dm_count,]
        head(no_dm_avgs)
        dim(no_dm_avgs)

        # get the dm regions
        get_dm_regions <- function(dm_value, regions_per_dm_value){
            dm_data <- load_simulated_data(dm_value) 
            
            dm_avgs <- dm_data$avgs[1:regions_per_dm_value,]
            dm_rpcs <- dm_data$rpcs[1:regions_per_dm_value,]

            list(avgs=dm_avgs, rpcs=dm_rpcs)
        }
        # determine balues for dm
        dm_range <- seq(0.00001, 0.1, 0.01)
        #dm_range <- seq(0.00001, 0.006, 0.0006)
        dm_range

        # compute the number of regions to take from each dm value
        (regions_per_dm_value <- dm_count / length(dm_range))

        dm_regions <- lapply(dm_range, get_dm_regions, regions_per_dm_value)
        dm_avgs <- do.call(rbind, lapply(dm_regions, function(x) x$avgs))
        dm_rpcs <- do.call(rbind, lapply(dm_regions, function(x) x$rpcs))
        dim(dm_avgs)
        dim(dm_rpcs)

        # combine no dm with dm data
        avgs <- rbind(no_dm_avgs, dm_avgs)
        rpcs <- rbind(no_dm_rpcs, dm_rpcs)
        dim(rpcs)
        dim(avgs)

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

        names(avgs_res)
        res <- avgs_res
        get_sig_counts <- function(res, st){
            names(res)
            results <- res$results
            head(results)
            tail(results)
            dim(results)

            # add the true label to the results
            results_labelled <- results %>%
                rownames_to_column('region') %>%
                mutate(pc_label = str_match(region, "-PC[0-9]*"),
                       pc_label = str_remove(pc_label, "-")) %>%
                mutate(region = str_remove(region, "-PC[0-9]*")) %>%
                mutate(percent_dm_difference = str_match(region, "_.*$"),
                       percent_dm_difference = str_remove(percent_dm_difference, "_"),
                       percent_dm_difference = as.numeric(percent_dm_difference)
                       ) %>%
                mutate(sig_dm = bf_pval < 0.05)
            head(results_labelled)
            table(results_labelled$percent_dm_difference)

            formatted_results <- results_labelled %>%
                select(region, pc_label, percent_dm_difference, sig_dm) %>%
                mutate(summary_type = st)
            head(formatted_results)

            formatted_results
        }

        avgs_labelled <- get_sig_counts(avgs_res, 'avgs')
        rpcs_labelled <- get_sig_counts(rpcs_res, 'rpcs')

        head(avgs_labelled)
        head(rpcs_labelled)
        table(rpcs_labelled$pc_label)

        true_labelled <- avgs_labelled %>%
            mutate(summary_type = 'True') %>%
            mutate(sig_dm = percent_dm_difference != 0)
        head(true_labelled)

        results_labelled <- rbind(avgs_labelled, rpcs_labelled) %>%
            rbind(true_labelled) %>%
            mutate(sig_label = ifelse(sig_dm, 'DM detected', 'No DM detected')) %>%
            mutate(percent_dm_difference = ifelse(percent_dm_difference > 1e-05,
                                                  percent_dm_difference - 0.00001,
                                                  percent_dm_difference)
                       ) %>%
            mutate(summary_type = ifelse(summary_type == 'rpcs', 'pcs', summary_type)) 
        head(results_labelled)
        table(results_labelled$summary_type)

        table(results_labelled$summary_type, results_labelled$sig_label)

        counts <- results_labelled %>%
            count(percent_dm_difference, summary_type, sig_label) %>%
            spread(summary_type, n, fill=0)
        head(counts)

        # compute the counts
        detailed_counts <- results_labelled %>%
            count(percent_dm_difference, sig_dm, summary_type)  %>%
            spread(summary_type, n, fill=0) %>%
            mutate(avg_err = True - avgs,
                   pc_err = True - pcs
            ) %>%
            filter(percent_dm_difference == 0 | sig_dm) %>%
            filter(percent_dm_difference > 1e-05) %>%
            mutate(avg_prop = avgs/True,
                   rpc_prop = pcs/True
            ) %>%
            select(percent_dm_difference, avg_prop, rpc_prop) %>%
            gather('summary_type', 'prop', -percent_dm_difference) %>%
            group_by(summary_type) %>%
            summarize(min=min(prop),
                      med=median(prop),
                      mean=mean(prop),
                      max=max(prop)
            )
        head(detailed_counts)



        summary_type_colors
        custom_st_colors <- summary_type_colors
        custom_st_colors['True'] <- 'grey'
        custom_st_colors

        head(results_labelled)
        unique(results_labelled$summary_type)

        p <- results_labelled %>%
            left_join(summary_type_map) %>%
            mutate(formatted_st = ifelse(!summary_type %in% c('avgs', 'pcs'), 
                                         as.character('True'), 
                                         as.character(formatted_st))) %>%
            ggplot() +
            geom_bar(aes(x=as.factor(percent_dm_difference), 
                         fill=formatted_st),
                     position=position_dodge2(preserve='single')
                     ) +
            facet_wrap(sig_label ~ ., nrow=2, scales='free_y') +
            theme_bw() +
            theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1),
                  strip.background=element_blank(), # remove facet strip background
                  panel.border=element_rect(color='black',
                                            fill=NA,
                                            linewidth=1
                                            )
                      ) +
            scale_fill_manual(values=custom_st_colors, name="Summary type") +
            xlab(paste("Percent methylation difference")) +
            ylab(paste("DM regions (BF p-value < 0.05)"))
        (savefile <- paste0(savedir, "mixed_layered_dm_results.png"))
        #(savefile <- paste0(savedir, "verySmall_mixed_layered_dm_results.png"))
        ggsave(p, file=savefile)

    }
}

run_mixed_dmp <- function(){
    # combine DM runs with no-DM runs to test sensitivity/specificity

    # select parameters
    percent_meth_low <- 0
    percent_meth_high <- 0.75

    num_sites <- 50
    num_samples <- 500
    percent_meth_difference <- 0.5

    #percent_sites_dm <- percent_meth_low
    load_simulated_data <- function(percent_sites_dm){
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
        
        names(results)

        # get labels for regions about DM
        simulations <- results$raw_data
        length(simulations)

        i <- 1
        check_region_dm_status <- function(i){
            simi <- simulations[[i]]
            if (all(is.na(simi))) return(NA)

            simi[[2]]

            # get the DM status of sites in this region
            dm_sites <- simi[[1]][,1:2]
            head(dm_sites)

            # count how many sites were DM
            tibble(num_dm_sites=sum(dm_sites$dmsites != 0),
                 num_sites_total=nrow(dm_sites))

        }
        dm_sites_list <- lapply(1:length(simulations), check_region_dm_status)
        dm_sites <- do.call(rbind, dm_sites_list) %>%
            mutate(region = row_number()) 
        head(dm_sites)

        summary(dm_sites$num_dm_sites)


        # grab the averages and rpcs
        avgs <- results$avgs
        rpcs <- results$rpcs

        rownames(avgs) <- paste(rownames(avgs), percent_sites_dm, sep='_')
        rownames(rpcs) <- paste(rownames(rpcs), percent_sites_dm, sep='_')

        list(avgs=avgs, rpcs=rpcs, dm_sites=dm_sites)
    }
    
    process_run <- function(){
        dmr_length <- num_sites / 2
        N <- 1000

        # get simulated data 
        simulated_data_low <- load_simulated_data(percent_meth_low)
        simulated_data_high <- load_simulated_data(percent_meth_high)

        # generate a mixed data frame
        n_low <- 400
        n_high <- 600

        # get the data from the simulations; mixing from the parameter runs
        avgs <- rbind(simulated_data_low$avgs[1:n_low,],
                          simulated_data_high$avgs[1:n_high,]
        )
        rpcs <- rbind(simulated_data_low$rpcs[1:n_low,],
                          simulated_data_high$rpcs[1:n_high,]
        )
        dm_sites <- rbind(simulated_data_low$dm_site[1:n_low,],
                          simulated_data_high$dm_sites[1:n_high,]
        ) %>%
            mutate(region = c(rownames(avgs)))
        head(dm_sites)


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

        names(avgs_res)
        res <- avgs_res
        get_sig_counts <- function(res, st){
            names(res)
            results <- res$results
            head(results)
            tail(results)
            dim(results)

            head(dm_sites)

            # add the true label to the results
            results_labelled <- results %>%
                rownames_to_column('region') %>%
                mutate(pc_label = str_match(region, "-PC[0-9]*"),
                       pc_label = str_remove(pc_label, "-")) %>%
                mutate(region = str_remove(region, "-PC[0-9]*")) %>%
                left_join(dm_sites) %>%
                mutate(true_dm = num_dm_sites > 0) %>%
                mutate(sig_dm = bf_pval < 0.05)
            head(results_labelled)

            formatted_results <- results_labelled %>%
                select(region, num_dm_sites, num_sites_total, true_dm, sig_dm) %>%
                mutate(summary_type = st)
            head(formatted_results)

            formatted_results
        }

        avgs_labelled <- get_sig_counts(avgs_res, 'avgs')
        rpcs_labelled <- get_sig_counts(rpcs_res, 'rpcs')

        head(avgs_labelled)
        head(rpcs_labelled)

        #results_labelled <- rbind(avgs_labelled, rpcs_labelled) %>%
            #mutate(outcome = ifelse(true_dm & sig_dm, 'tp', NA),
                   #outcome = ifelse(true_dm & !sig_dm, 'fn', outcome),
                   #outcome = ifelse(!true_dm & !sig_dm, 'tn', outcome),
                   #outcome = ifelse(!true_dm & sig_dm, 'fp', outcome),
            #)
        results_labelled <- rbind(avgs_labelled, rpcs_labelled) %>%
            separate(region, c('region', 'percent_sites_dm'), sep='_') %>%
            mutate(true_dm = percent_sites_dm == 0.75) %>%
            mutate(outcome = ifelse(true_dm & sig_dm, 'tp', NA),
                   outcome = ifelse(true_dm & !sig_dm, 'fn', outcome),
                   outcome = ifelse(!true_dm & !sig_dm, 'tn', outcome),
                   outcome = ifelse(!true_dm & sig_dm, 'fp', outcome),
            )
        head(results_labelled)

        table(results_labelled$outcome, results_labelled$summary_type)

        p <- ggplot(aes())

    }
}

main <- function(runnum, N){
    randomized_pheno <- TRUE
    
    # run DM analysis
    run <- runs[212,] # for testing interactively
    process_run <- function(run){
        num_sites <- run[['num_sites']] %>% as.numeric()
        num_samples <- run[['num_samples']] %>% as.numeric()
        percent_meth_difference <- run[['percent_meth_difference']] %>% as.numeric()
        percent_sites_dm <- run[['percent_sites_dm']] %>% as.numeric()

        dmr_length <- num_sites / 2
        N <- 1000

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
        
        names(results)

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

        if (randomized_pheno){
            # we want to randomize phenotype to test 
            # number of false positives 
            set.seed(71446174)
            rand_pheno <- pheno
            rand_pheno$status <- sample(pheno$status, nrow(pheno), replace=FALSE)
            head(rand_pheno)
            table(pheno$status, rand_pheno$status)
            pheno <- rand_pheno
        }
        
        avgs_res <- run_dmp(avgs, pheno)
        rpcs_res <- run_dmp(rpcs, pheno)


        dm_res <- list(avgs_res=avgs_res, rpcs_res=rpcs_res)
        
        # compare methods
        head(avgs_res)
        names(avgs_res)
        
        #dm_res <- avgs_res
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


        # store the number of sig results for avgs and rpcs
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

        names(dm_res)
        test <- dm_res$num_sig
        class(test)
        names(test)
        head(test)

        (filename <- paste0("DM_results",
                            "_numSites", num_sites,
                            "_numSamples", num_samples,
                            "_pct_meth_diff", percent_meth_difference,
                            "_dmr_length", dmr_length,
                            "_pct_sites_dm", percent_sites_dm,
                            "_N", N))
        if (randomized_pheno){
            (filename <- paste0("randomizedPheno_DM_results",
                                "_numSites", num_sites,
                                "_numSamples", num_samples,
                                "_pct_meth_diff", percent_meth_difference,
                                "_dmr_length", dmr_length,
                                "_pct_sites_dm", percent_sites_dm,
                                "_N", N))
        }
        (savefile <- paste0(savedir, filename, ".rds"))
        saveRDS(dm_res, savefile)
    }
    
    
    # round 1
    num_sites_range <- c(20,50)
    #num_samples_range <- c(50,500,5000)
    num_samples_range <- c(50,500)
    #num_samples_range <- c(5000)
    #percent_meth_difference_range <- seq(0.1, 1, 0.1)
    percent_meth_difference_range <- seq(0.00001, 0.1, 0.01)
    #percent_meth_difference_range <- seq(0.00001, 0.006, 0.0006)
    #percent_meth_difference_range <- seq(0.01, 0.1, 0.01)
    percent_sites_dm_range <- seq(0.25,1,0.25)

    # round 2
    #num_sites_range <- c(20,50)
    #num_samples_range <- c(50,500)
    #percent_meth_difference_range <- seq(0.01, 0.2, 0.01)
    #percent_sites_dm_range <- seq(0,0.75,0.25)

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
