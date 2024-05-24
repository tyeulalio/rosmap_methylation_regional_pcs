library(tidyverse)

# check rosmap data for ranges of simulation parameters

savedir <- paste0("../../../output/parameters/")
dir.create(savedir)

check_percent_meth_difference <- function(){
    # check rosmap for different in methylation between groups
    # load full gene cpg differential methylation results
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/04_dmp_analysis/01_dmp_results/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/global_pcs/protect_global_pcs/nice_results/full_gene_full_DM_results.csv")
    dm_results <- read_csv(datafile)
    head(dm_results)

    # filter for cpgs only
    unique(dm_results$trait)
    ct <- 'astro'
    tr <- 'ceradsc'
    cpg_results <- dm_results %>%
        filter(summary_type == 'cpgs',
               trait == tr,
               cell_type == ct
        )
    head(cpg_results)
    dim(cpg_results)
    table(cpg_results$bh_pval < 0.05, cpg_results$cell_type)


    if (ct != 'bulk'){
        # load the raw beta values for astrocytes, cpgs
        datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/02_deconvolution/01_deconvolved_data/restimate_proportions/c2_covariates/tca_tensor.rds")
        tca_tensor <- readRDS(datafile)

        # get ordered cell types
        datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/02_deconvolution/01_deconvolved_data/restimate_proportions/c2_covariates/tca_out.rds")
        tca_out <- readRDS(datafile)

        names(tca_tensor) <- colnames(tca_out$W)

        betas_matrix <- tca_tensor[[ct]]
        head(betas_matrix[,1:10])

        print("scaling betas after TCA")
        betas_matrix[betas_matrix < 0] <- 0
        betas_matrix[betas_matrix > 1] <- 1
        betas_matrix <- (betas_matrix * (ncol(betas_matrix) - 1) + 0.5) / ncol(betas_matrix)
        head(betas_matrix)
        betas <- betas_matrix
    }
    if (ct == 'bulk'){
        datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/02_deconvolution/01_deconvolved_data/restimate_proportions/c2_covariates/formatted_betas.rds")
        betas <- readRDS(datafile)
        head(betas)
    }

    # load the formatted traits
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/02_deconvolution/02_svd_batch_correction/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/astro_formatted_clinical.rds")
    traits <- readRDS(datafile)
    head(traits)


    # need a mapping from hg19 to hg38 cpg names
    head(betas)
    head(cpg_results)

    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/03_summarised_genes/01_summarised_genes/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/", ct, "_cpg_position_map.rds")
    cpg_map <- readRDS(datafile)
    head(cpg_map)

    hg38_betas_cpgs <- data.frame(betas_hg19_cpgs = rownames(betas)) %>%
        mutate(cpg_gr37 = betas_hg19_cpgs) %>%
        left_join(cpg_map[c('cpg_gr37', 'cpg')]) %>%
        rename(cpg_gr38=cpg) 
    head(hg38_betas_cpgs)
    dim(hg38_betas_cpgs)

    length(intersect(hg38_betas_cpgs$betas_hg19_cpgs, rownames(betas)))
    length(intersect(hg38_betas_cpgs$cpg_gr38, cpg_results$feature_id))

    # figure out the samples that we're using
    # compare cerad group 1 to 4 for largest differences
    if (tr == 'ceradsc'){
        traits_filtered <- traits %>%
            filter(ceradsc %in% c(1,4)) %>%
            select(ceradsc) %>%
            arrange(ceradsc) 
    }
    if (tr == 'new_diag'){
        head(traits)
        unique(traits$new_diag)
        traits_filtered <- traits %>%
            filter(new_diag %in% c('AD', 'Control')) %>%
            select(new_diag) %>%
            arrange(new_diag)
    }
    traits_filtered <- traits_filtered %>%
        mutate(individual_id = rownames(traits_filtered))
    head(traits_filtered)
    table(traits_filtered$ceradsc)
    table(traits_filtered$new_diag)

    # filter the betas, get average difference between group 1 and 4
    head(betas)
    betas_filt <- betas[,rownames(traits_filtered)]
    dim(betas_filt)
    dim(traits_filtered)
    head(betas_filt)

    betas_long <- betas_filt %>%
        as.data.frame() %>%
        rownames_to_column('cpg_gr37') %>%
        gather('individual_id', 'meth', -cpg_gr37)
    head(betas_long)

    head(traits_filtered)

    # attach the trait
    betas_trait <- betas_long %>%
        left_join(traits_filtered)
    head(betas_trait)
    summary(betas_trait$meth)


    if (tr == 'ceradsc'){
        betas_diff <- betas_trait %>%
            group_by(ceradsc, cpg_gr37) %>%
            summarize(mean_meth = mean(meth))
        wide_diff <- betas_diff %>%
            mutate(ceradsc = paste0("ceradsc", ceradsc)) %>%
            spread(ceradsc, mean_meth) %>%
            mutate(beta_diff = ceradsc4 - ceradsc1)
    } 
    if (tr == 'new_diag'){
        betas_diff <- betas_trait %>%
            #group_by(ceradsc, cpg_gr37) %>%
            group_by(new_diag, cpg_gr37) %>%
            summarize(mean_meth = mean(meth))
        wide_diff <- betas_diff %>%
            spread(new_diag, mean_meth) %>%
            mutate(beta_diff = AD - Control)
    }
    head(wide_diff)

    # mark whever these were found s DM or not
    # add hg38 names
    head(cpg_results)
    sub_dm <- cpg_results %>%
        select(cpg_gr38=feature_id, bh_pval) %>%
        mutate(dm = bh_pval < 0.05) %>%
        left_join(hg38_betas_cpgs) 
    head(sub_dm)

    # join the dm with the betas diffs
    diffs_dm <- wide_diff %>%
        inner_join(sub_dm)
    head(diffs_dm)
    dim(diffs_dm)

    table(diffs_dm$dm)
    summary(abs(diffs_dm$beta_diff))
    summary(abs(diffs_dm[diffs_dm$dm,]$beta_diff))
    summary(abs(diffs_dm[!diffs_dm$dm,]$beta_diff))


    p <- diffs_dm %>%
        ggplot() +
        geom_histogram(aes(x=abs(beta_diff), fill=dm),
                       bins=100
                       #alpha=0.5
        ) +
        facet_wrap(. ~ dm, nrow=2, scales='free_y') +
        coord_cartesian(xlim=c(0,0.04)) +
        theme_bw()
    (savefile <- paste0(savedir, tr, "_", ct, "_density_plot_dm.png"))
    ggsave(p, file=savefile)

    # what is the proportion DM vs not DM at each quantile?
    head(diffs_dm)
    quants <- quantile(abs(diffs_dm$beta_diff), seq(0, 1, length.out=30))
    quants

    quants_dm <- diffs_dm %>%
        mutate(quant_groups = cut(abs(beta_diff), breaks=quants, include.lowest =TRUE))
    head(quants_dm)

    quant_counts <- quants_dm %>%
        count(dm, quant_groups) %>%
        mutate(dm = ifelse(dm, 'dm', 'not_dm')) %>%
        spread(dm, n, fill=0) %>%
        mutate(prop_dm = dm / (dm+not_dm)) %>%
        mutate(quant_order = as.numeric(quant_groups))
    head(quant_counts)

    p <- quant_counts %>%
        ggplot() +
        geom_bar(aes(x=fct_reorder(quant_groups, quant_order), y=prop_dm),
                 stat='identity'
        ) + 
        #geom_hline(yintercept=0.05, color='red') +
        theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))
    (savefile <- paste0(savedir, tr, "_proportion_plot.png"))
    ggsave(p, file=savefile, width=12, height=5)

}

main <- function(){
    # check percent methylation difference
    check_percent_meth_difference()
}
