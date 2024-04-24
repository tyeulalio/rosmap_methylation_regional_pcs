library(ggpubr)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

# plot the disease-relevant of PC, differential methylation figure


savedir <- paste0("/path/to/output/")
dir.create(savedir)

combine_plots <- function(plot1, plot2, plot3){
    
    p <- plot1
    format_plot <- function(p){
        p + theme(text=element_text(size=20),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  legend.text=element_text(size=18),
                  legend.title=element_text(size=20),
                  plot.margin=margin(t=10,
                                     r=20,
                                     b=10,
                                     l=10
                  )
        ) 
    }

    format_plot1 <- function(p){
        tmp <- summary_type_colors
        names(tmp) <- c('Avgs', 'rPCs', 'CpGs', 'white')
        p + theme(legend.position=c(0.65,0.80)) +
            scale_fill_manual(values=tmp,
                              name="Summary type"
            )
    }
    format_plot2 <- function(p){
        p + guides(linetype=guide_legend(title=paste0("Summary type")))
    }
    format_plot3 <- function(p){
        tmp <- summary_type_colors
        tmp[4] <- 'darkgrey'
        names(tmp) <- c('Avgs only', 'rPCs only', 'CpGs', 'Both')
        p + theme(
                  axis.text.x=element_text(size=20),
                  axis.text.y=element_text(size=24)
        ) +
            scale_fill_manual(values=tmp, 
                              name="Differential\nmethylation",
                              na.value='white')
    }

    formatted_p1 <- format_plot(plot1)
    formatted_p1 <- format_plot1(formatted_p1)

    formatted_p2 <- format_plot(plot2)
    formatted_p2 <- format_plot2(formatted_p2)

    formatted_p3 <- format_plot(plot3)
    formatted_p3 <- format_plot3(formatted_p3)

    # create the top row of the plot
    top_p <- ggarrange(formatted_p1, 
                       formatted_p2, 
                       labels=c('a','b'),
                       nrow=1,
                       widths=c(0.75, 1)
    )

    # combine all plots
    combined_p <- ggarrange(top_p,
                            formatted_p3,
                            nrow=2,
                            labels=c('','c'),
                            heights=c(1,0.75)
    )

    (savefile <- paste0(savedir, "dm_results_figure.png"))
    ggsave(combined_p, file=savefile, width=13, height=9)
}

main <- function(){
    # load the upset plot 
    datadir <- paste0("/path/to/new_rosmap/output/04_dmp_analysis/01_dmp_results/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/global_pcs/protect_global_pcs/plots/bh_matched_array/")
    datafile <- paste0(datadir, "ceradsc_celltype_intersect_barplot.rds")
    plot1 <- readRDS(datafile)

    # load ad gene tile plot
    datafile <- paste0(datadir, "ceradsc_bh_quant0.95_opentargets_focusgenes_tileplot.rds")
    plot3 <- readRDS(datafile)

    # load interaction model plot
    datadir <- paste0("/path/to/new_rosmap/output/04_dmp_analysis/02_processed_dmp_results/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/plots/bh_matched_array/")
    datafile <- paste0(datadir, "ceradsc_full_gene_pval_adScore_scatter.rds")
    plot2 <- readRDS(datafile)

    # combine pltos
    combine_plots(plot1, plot2, plot3)
}
