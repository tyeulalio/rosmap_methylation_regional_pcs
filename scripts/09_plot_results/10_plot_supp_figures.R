library(ggplot2)
library(ggpubr)
#library(gridExtra)
library(tidyverse)

# plot supplementary figures that need formatting
savedir <- paste0("/path/to/output/")
dir.create(savedir)

plot_supp4 <- function(){
    # load plots
    datafile <- paste0("/path/to/new_rosmap/output/04_dmp_analysis/01_dmp_results/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/global_pcs/protect_global_pcs/plots/bh_matched_array/ceradsc_bh_fullGene_dm_genes_counts_barplot.rds")
    gene_plot <- readRDS(datafile)

    datafile <- paste0("/path/to/new_rosmap/output/04_dmp_analysis/01_dmp_results/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/global_pcs/protect_global_pcs/plots/bh_matched_array/ceradsc_feature_proportions_comparison.rds")
    feat_plot <- readRDS(datafile)

    # format the plot
    feat_plot2 <- feat_plot +
        theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1),
              legend.position='top',
              legend.spacing.x=unit(0.75, 'cm')
        )
    gene_plot2 <- gene_plot +
        theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1)
        )

    # combine plots
    p <- ggarrange(
                    feat_plot2, 
                   gene_plot2, 
                       labels=c('A','B'),
                       nrow=1,
                       common.legend=TRUE,
                       widths=c(0.75, 1)
    ) 
    (savefile <- paste0(savedir, "supp_fig4.png"))
    ggsave(p, file=savefile, width=10, height=7)

}

main <- function(){
    # plot supp figure 4
    plot_supp4()
}
