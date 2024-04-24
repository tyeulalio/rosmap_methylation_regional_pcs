library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

# plot manuscript figure with cell type and summary type sections

savedir <- paste0("/path/to/output/")
dir.create(savedir)


load_celltype_props <- function(){
    datafile <- paste("/path/to/new_rosmap/output/02_deconvolution/01_deconvolved_data/restimate_proportions/c2_covariates/include_study/plots/episcore_celltype_proportions_boxplot.rds")
    celltype_props_fig <- readRDS(datafile)
    celltype_props_fig
}

load_celltype_umap <- function(){
    datafile <- paste0("/path/to/new_rosmap/output/02_deconvolution/03_celltype_umap/restimate_proportions/c2_covariates/mvals/remove_age/noCleaningSex/cleaned_celltype_umap_plot.rds")
    celltype_umap_fig <- readRDS(datafile)
    celltype_umap_fig
}

combine_plots <- function(celltype_props_fig, celltype_umap_fig, summarytype_fig, cpg_pc_fig){ 

    # format the plots
    common_theme <- theme(text=element_text(size=20),
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                          strip.background=element_blank(),
                          panel.border=element_rect(color='black', 
                                                    fill=NA,
                                                    linewidth=1
                          )
    ) 

    a_formatted <- celltype_props_fig +
        common_theme +
        theme(axis.text.x = element_text(angle=25, hjust=1, vjust=1)) +
        ylab("Estimated proportions")
    b_formatted <- celltype_umap_fig +
        common_theme +
        guides(color = guide_legend(override.aes=list(size=10))) +
        theme(legend.spacing.x = unit(0.4, 'cm'))
    #c_formatted <- summarytype_fig +
        #common_theme +
        #scale_color_manual(values=summary_type_colors)
    d_formatted <- cpg_pc_fig +
        common_theme +
        scale_color_manual(values=cell_type_colors) 

    # combine the plots into one
    p <- ggarrange(a_formatted, 
                   b_formatted,
                   #c_formatted,
                   d_formatted,
                   #labels = c("A", "B", "C", "D"),
                   labels = 'auto',
                   ncol=3,
                   nrow=1,
                   common.legend = TRUE,
                   vjust=1,
                   font.label=list(size=20)
    ) 
    (savefile <- paste0(savedir, "celltype_summarytype_figure.png"))
    ggsave(savefile, p, height=5, width=13)

}

load_summarized_genes_data <- function(cell_type, region_type){
    region_type <- 'full_gene'
    #cell_type <- 'astro'
    cleand_str = "cleaned_global_"

    print(paste("Loading", cell_type, region_type))

    # load the summary type data
    datadir <- paste0("/path/to/new_rosmap/output/03_summarised_genes/01_summarised_genes/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/")

    # reading in data

    # read PCs
    datafile <- paste0(datadir, cleand_str, region_type, "_", cell_type,
                      "_pcs.rds")
    datafile
    if (region_type != 'full_gene'){
    pcs <- readRDS(datafile) %>%
      select(-num_cpgs, -gv_dim, -mp_dim, -error, -duration_secs, -percent_variance)
    } else{
    pcs <- readRDS(datafile) %>%
      select(-percent_variance_explained)
    }

    head(pcs)

    # read averages
    datafile <- paste0(datadir, cleand_str, region_type, "_", cell_type,
                     "_avgs.rds")
    avgs <- readRDS(datafile) %>%
        select(-duration_secs)
    head(avgs)

    # read in cpgs - skip this
    #datafile <- paste0(datadir, cleand_str,
                     #cell_type, "_lifted_hg38_mvals.rds"
    #)
    #cpgs <- readRDS(datafile)
    #head(cpgs)

    # read in position map
    #datafile <- paste0(datadir, cell_type, "_cpg_position_map.rds")
    #pos_map <- readRDS(datafile) %>%
        #select(group_name=seqnames, cpg)
    #head(pos_map)

  # replace probes with positions
  #head(cpgs)
  #formatted_cpgs <- cpgs %>%
    #as.data.frame()
  #head(formatted_cpgs)


  # read in cpg map
  #datafile <- paste0(datadir, region_type, "_",
                     #cell_type, "_gene_cpg_map.rds")
  #cpg_map <- readRDS(datafile)
  #head(cpg_map)


    #dat <- list(pcs=pcs, cpgs=formatted_cpgs, avgs=avgs, cpg_map=cpg_map)
    dat <- list(pcs=pcs, avgs=avgs)
    return(dat)
}

# get genes and their symbols
get_st_genes <- function(matched_dat){
    genes <- data.frame(gene_id=rownames(matched_dat$avgs)) %>%
        separate(gene_id, c('gene', 'gene_vers'), remove=FALSE, sep='\\.')
    head(genes)

    # load gene symbols
    datafile <- paste0("/path/to/data/full_biomart_table.csv")
    gene_symbols <- read_csv(datafile) %>%
        rename(gene_id=Ensembl_ID,
               symbol=Symbol
        ) %>%
        separate(gene_id, c('gene', 'gene_vers'), remove=FALSE, sep='\\.')
    head(gene_symbols)

    # load ad gene annotations
    datafile <- paste0("/path/to/data/gene_sets/AD_opentargets.csv")
    ad_genes <- read_csv(datafile)
    head(ad_genes)

    # attach the gene 
    gene_annots <- genes %>%
        left_join(gene_symbols) %>%
        left_join(ad_genes[c('symbol', 'overallAssociationScore')])
    head(gene_annots)

    gene_annots
}

make_summarytype_fig <- function(){
    # make summary type figure

    # load the data for this figure
    cell_type = 'astro'
    region_type = 'full_gene'
    matched_dat <- load_summarized_genes_data(cell_type, region_type)

    # load some phenotype data
    datafile <- paste0("/path/to/new_rosmap/output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds")
    clinical <- readRDS(datafile)
    head(clinical)

    # get genes and their symbols
    genes <- get_st_genes(matched_dat)
    head(genes)

    # load differential methylation results
    datafile <- paste0("/path/to/new_rosmap/output/04_dmp_analysis/01_dmp_results/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/global_pcs/protect_global_pcs/combined_dm_results.rds")
    dm_results <- readRDS(datafile)
    head(dm_results)
    unique(dm_results$feature_id)

    # find gene DM in pcs but not averages
    interest_dm <- dm_results %>%
        #mutate(sig = bh_pval < 0.05) %>%
        #separate(feature_id, c('gene_id', 'pc'), sep='-', fill='right') %>%
        rename(pc=feature_id) %>%
        select(summary_type, cell_type, gene_id, bh_pval) %>%
        group_by(summary_type, gene_id, cell_type) %>%
        #summarize(sig = sum(sig) > 0) %>%
        summarize(pval = min(bh_pval)) %>%
        filter(!is.na(summary_type)) %>%
        spread(summary_type, pval) %>%
        #filter(pcs, !avgs)
        filter(pcs < 0.05, avgs > 0.05) %>%
        mutate(diff = avgs - pcs) %>%
        arrange(-diff) %>%
        ungroup() 
    head(interest_dm)
    dim(interest_dm)

    # mark whether genes are DM 
    head(genes)
    genes2 <- genes %>%
        mutate(dm_gene = gene %in% interest_dm$gene_id) %>%
        select(-gene_id) %>%
        rename(gene_id=gene) %>%
        inner_join(interest_dm) %>%
        arrange(-overallAssociationScore, -diff) %>%
        filter(!is.na(cell_type)) 
    head(genes2)
    unique(genes2$symbol)[1:20]

    # select the gene that we're interested in
    gene_symbol <- "PSENEN"

    genes <- genes2
    genes %>%
        filter(symbol == gene_symbol)

    ct <- 'oligo_opc'


    filtered_genes <- genes %>%
        filter(symbol == gene_symbol,
               cell_type == ct
        )
    head(filtered_genes)
    (gene_id <- filtered_genes[['gene_id']])

    gid = gene_id
     dm_results %>%
        mutate(sig = bh_pval < 0.05) %>%
        #separate(feature, c('gene_id', 'pc'), sep='-', fill='right') %>%
        select(summary_type, cell_type, sig, gene_id, pc) %>%
        filter(gene_id == gid
               #cell_type == ct
        )


    matched_dat <- load_summarized_genes_data(ct, region_type)

    # grab the gene averages
    avg <- matched_dat$avgs %>%
        rownames_to_column('gene') %>%
        #filter(gene == gene_id) %>%
        filter(grepl(gene_id, gene)) %>%
        gather(key='sample_name', value='avg', -gene) %>%
        select(sample_name, avg) %>%
        column_to_rownames('sample_name')
    head(avg)

    # scale averages
    avg <- scale(avg)
    avg <- avg %>%
        as.data.frame() %>%
        rownames_to_column('sample_name')
    head(avg)
  
  
    # grab gene apc
    gene_id
    head(matched_dat$pcs)
    pcs <- matched_dat$pcs %>%
        rownames_to_column('gene') %>%
        separate('gene', into=c('gene', 'pc'), sep='-', fill='right') %>%
        #filter(gene == gene_id) %>%
        filter(grepl(gene_id, gene)) %>%
        gather(key='sample_name', value='apc', -gene, -pc) %>%
        spread(key='pc', value='apc') %>%
        select(-gene) %>%
        column_to_rownames('sample_name')
    head(pcs)
  
    pcs <- scale(pcs)
    pcs <- pcs %>%
        as.data.frame() %>%
        rownames_to_column('sample_name')
    head(pcs)
  
    # combine the data
    head(pcs)
    head(avg)
  
    head(clinical)
    sub_clinical <- clinical %>%
        select(sample_name=wgs.individualID, ceradsc)
    head(sub_clinical)
  
    combined_meth <- inner_join(pcs, avg) %>%
        left_join(sub_clinical)
    head(combined_meth)

    # get pvalue labels
    head(genes2)
    labels <- genes2 %>%
        filter(symbol == gene_symbol, 
               cell_type == ct
        ) %>%
        select(Avgs=avgs, rPCs=pcs) %>%
        gather('formatted_st', 'pval')
    head(labels)
  
    gene_symbol
    p <- combined_meth %>%
        gather('summary_type', 'meth', PC1, avg) %>%
        mutate(summary_type = ifelse(summary_type == 'PC1', 'pcs', 'avgs')) %>%
        left_join(summary_type_map) %>%
        ggplot(aes(x=as.factor(ceradsc), y=meth)) +
        # geom_point(aes(x=PC1, y=avg, color=as.factor(ceradsc))) +
        # geom_point() +
        geom_violin(linewidth=1, fill='lightgrey', color='lightgrey') +
        geom_boxplot(linewidth=1, color='darkgrey', width=0.5) +
        geom_jitter(aes(color=formatted_st), alpha=0.8) +
        geom_label(data=labels, aes(x=3.75, y=2.75, label=paste0("BH p = ", round(pval,3)))) +
        facet_wrap(formatted_st ~ .) +
        geom_smooth(aes(x=ceradsc, y=meth),
                    color='blue', 
                    method='lm', se=FALSE) +
        scale_color_manual(values=summary_type_colors, 
                           name="Summary type"
                           #breaks=c("Avgs", "PCs")
                           ) +
        coord_cartesian(ylim=c(-3,3)) +
        theme_bw() +
        theme(text = element_text(size=20),
              legend.position='none',
              legend.justification='right',
              strip.background=element_rect(fill='white', color='white')) +
        xlab(paste("CERAD score")) +
        ylab(expression(paste(italic("PSENEN"), " normalized methylation")))

    (savefile <- paste0(savedir, gene_symbol, "_summarytype_violin.png"))
    ggsave(p, file=savefile)
  
    (savefile <- paste0(savedir, gene_symbol, "_summarytype_violin.rds"))
    saveRDS(p, savefile)

    return(p)
}

make_cpg_pc_fig <- function(){
    # make the cpg-pc count scatter plot
    datadir <- paste0("/path/to/new_rosmap/output/03_summarised_genes/02_check_genes/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/")

    # read in the counts
    list.files(datadir)
    datafile <- paste0(datadir, "pc_cpg_counts.rds")
    counts <- readRDS(datafile)
    head(counts)

    # plot cpgs with pcs for full gene, bulk
    fitlm <- lm(pc_count ~ cpg_count, data=counts)

    formatted_counts <- counts %>%
        left_join(cell_type_map) %>%
        filter(region_type == 'full_gene') 
    head(formatted_counts)
    unique(formatted_counts$cell_type)

    p <- ggplot(formatted_counts,
                aes(x=cpg_count, y=pc_count, color=formatted_ct)) +
        geom_point(
                   size=3,
                   show.legend=FALSE
        ) +
        geom_smooth(method='loess', se=FALSE) +
        theme_bw() +
        xlab(paste("CpGs per gene region")) +
        ylab(paste("rPCs per gene region")) +
        scale_color_manual(values=cell_type_colors) +
        guides(color=guide_legend(title=paste("Cell type")))
    return(p)
}

main <- function(){
    # load cell type proportions figure
    celltype_props_fig <- load_celltype_props()

    # load cell type umap figure
    celltype_umap_fig <- load_celltype_umap()

    # plot summary types
    summarytype_fig <- make_summarytype_fig()

    # plot cpg-pc plot
    cpg_pc_fig <- make_cpg_pc_fig()

    # combine plots
    combine_plots(celltype_props_fig, celltype_umap_fig, summarytype_fig, cpg_pc_fig)
}
