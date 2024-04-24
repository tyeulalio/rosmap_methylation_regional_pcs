# get counts across QTL analysis
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)


# creates the supplementary table with counts for GWAS integration results

savedir <- paste0("/path/to/ouput/")
dir.create(savedir)

plotdir <- paste0(savedir, "02_get_counts_plots/")
dir.create(plotdir)


get_total_counts <- function(runname){
    # get total counts of genes x variants tested
    datafile <- paste0("/path/to/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/01_torus_input/with_proportions/", runname, "_fastqtl_single_snp_output.tsv")

    # count number of gene x variants tested
    # this is actually feature-level for PCs (i.e. each PC for a gene counts as one)
    (cmd <- paste0("cut -d' ' -f1-2 ", datafile, " | uniq | wc -l"))
    #(cmd <- paste0("cut -d' ' -f1-2 ", datafile))
    (feature_var_count <- system(cmd, intern=TRUE))
    #feature_vars <- system(cmd, intern=TRUE) %>%
        #as.data.frame()
    #colnames(feature_vars) <- c('feature_var')
    #head(feature_vars)

    #formatted_tested <- feature_vars %>%
        #separate(feature_var, c('feature_id', 'variant_id'), sep=' ') %>%
        #separate(feature_id, c('gene_id', 'pc'), sep='-', fill='right')
    #head(formatted_tested)

    # count number of unique genes (feature-level)
    (cmd <- paste0("cut -d' ' -f1 ", datafile, " | uniq | wc -l"))
    (feature_count <- system(cmd, intern=TRUE))

    # count the number of genes for pcs
    (cmd <- paste0("cut -d' ' -f1 ", datafile))
    features <- system(cmd, intern=TRUE)
    head(features)

    genes <- str_remove(features, "-PC[0-9]*")
    head(genes)

    (gene_count <- length(unique(genes)))
    

    # count number of unique variants
    (cmd <- paste0("cut -d' ' -f2 ", datafile, " | uniq | wc -l"))
    (var_count <- system(cmd, intern=TRUE))

    # combine results
    data.frame(count_type=c('featureVar_0tested', 'feature_0tested', 'gene_0tested', 'var_0tested'),
               count = c(feature_var_count, feature_count, gene_count, var_count)
    )
}

get_sig_mapped_qtls <- function(old_runname){ 
    # get sig qtls after qtl mapping
    datafile <- paste0("/path/to/new_rosmap/output/05_qtl_analysis/12_map_qtls/with_proportions/",
                       old_runname, "/cis_qtl.signif_pairs.csv")
    sig_qtls <- read_csv(datafile) %>%
        separate(phenotype_id, c('phenotype_id', 'pc'), sep='-', fill='right')
    head(sig_qtls)

    # count gene x variant pairs
    gene_vars <- sig_qtls %>%
        select(phenotype_id, variant_id) %>%
        unique()
    head(gene_vars)
    (gene_var_count <- nrow(gene_vars))

    # get feature x variant pairs
    feature_vars <- sig_qtls %>%
        select(phenotype_id, pc, variant_id) %>%
        unique()
    head(feature_vars)
    (feature_var_count <- nrow(feature_vars))

    # count unique genes
    genes <- sig_qtls %>%
        select(phenotype_id) %>%
        unique()
    (gene_count <- nrow(genes))
    
    # count unique features
    features <- sig_qtls %>%
        select(phenotype_id, pc) %>%
        unique()
    (feature_count <- nrow(features))

    # count unique variants
    vars <- sig_qtls %>%
        select(variant_id) %>%
        unique()
    (var_count <- nrow(vars))

    
    data.frame(count_type=c('geneVar_1qtl', 'featureVar_1qtl', 'gene_1qtl', 'feature_1qtl', 'var_1qtl'),
               count = c(gene_var_count, feature_var_count, gene_count, feature_count, var_count)
    )
}

get_finemapped_qtls <- function(runname){
    # get counts of significant finemapped qtls
    datafile <- paste0("/path/to/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/05_fastqtl_qtl_annotations/with_proportions/", runname, "/fastenloc.qtl.annotation.vcf.gz")
    finemapped <- read_table(datafile, col_names=c('chrom', 'pos', 'variant_id', 'ref', 'alt', 'info'))
    head(finemapped)

    formatted_finemapped <- finemapped %>%
        separate_rows(info, sep='\\|') %>%
        separate(info, c('feature_id', 'cluster_info', 'num_snps'), sep=':') %>%
        separate(feature_id, c('gene_id', 'pc'), sep='-PC', fill='right') %>%
        mutate(num_snps = as.numeric(str_remove(num_snps, ']'))) %>%
        separate(cluster_info, c('cluster_info', 'cluster_pip'), sep='\\[') %>%
        mutate(cluster_pip = as.numeric(cluster_pip)) %>%
        filter(cluster_pip > 0.5)
    head(formatted_finemapped)

    # count gene x variant pairs
    gene_vars <- formatted_finemapped %>%
        select(gene_id, variant_id) %>%
        unique()
    head(gene_vars)
    (gene_var_count <- nrow(gene_vars))
    
    # count feature x variant pairs
    feature_vars <- formatted_finemapped %>%
        select(gene_id, pc, variant_id) %>%
        unique()
    head(feature_vars)
    (feature_var_count <- nrow(feature_vars))

    # count unique genes
    genes <- formatted_finemapped %>%
        select(gene_id) %>%
        unique()
    (gene_count <- nrow(genes))

    # count unique features
    features <- formatted_finemapped %>%
        select(gene_id, pc) %>%
        unique()
    (feature_count <- nrow(features))

    # count unique variants
    vars <- formatted_finemapped %>%
        select(variant_id) %>%
        unique()
    (var_count <- nrow(vars))

    head(formatted_finemapped)
    (cluster_size = summary(formatted_finemapped$num_snps))

    data.frame(count_type=c('geneVar_2finemapped', 'featureVar_2finemapped', 'gene_2finemapped', 
                            'feature_2finemapped', 'var_2finemapped',
                            'cluster_finemapped_min', 'cluster_finemapped_mean', 
                            'cluster_finemapped_median', 'cluster_finemapped_max'),
               count = c(gene_var_count, feature_var_count, gene_count,
                         feature_count, var_count,
                         cluster_size[[1]], cluster_size[[4]], 
                         cluster_size[[3]], cluster_size[[6]]
               )
    )
}

get_colocalization_counts <- function(runname, coloc_thresh){
    # get counts of colocalizations
    datafile <- paste0("/path/to/new_rosmap/output/06_fine_mapping_colocalization/04_fastenloc/01_fastenloc_output/with_proportions/", runname, "/wightman/fastqtl_ld.enloc.gene.out")
    gene_colocalizations <- read_table(datafile)
    head(gene_colocalizations)

    # filter for high colocalization probabilities
    sig_colocalizations <- gene_colocalizations %>%
        filter(GLCP > coloc_thresh) %>%
        separate(Gene, c('gene_id', 'pc'), sep='-PC', fill='right')
    sig_colocalizations
    (num_features <- nrow(sig_colocalizations))
    (num_genes <- unique(length((sig_colocalizations$gene_id))))

    (final_out <- data.frame(gene_3colocalized_GLCP = num_genes))

    # check other files for more info on variants per cluster and whatever
    datafile <- paste0("/path/to/new_rosmap/output/06_fine_mapping_colocalization/04_fastenloc/01_fastenloc_output/with_proportions/", runname, "/wightman/fastqtl_ld.enloc.sig.out")
    gene_var_colocalizations <- read_table(datafile)
    head(gene_var_colocalizations)

    # filter for genes with significant GLCPs
    sig_gene_var <- gene_var_colocalizations %>%
        separate(Signal, c('qtl_signal', 'gwas_signal'), sep="\\(@\\)") %>%
        separate(qtl_signal, c('feature_id', 'qtl_cluster_id'), sep=":") %>%
        separate(feature_id, c('gene_id', 'pc'), sep='-PC', remove=FALSE, fill='right') %>%
        separate(gwas_signal, c('gwas_ld_block', 'gwas_cluster_id'), sep=":") %>%
        filter(gene_id %in% sig_colocalizations$gene_id) %>%
        mutate(Num_SNP = as.numeric(Num_SNP))
    head(sig_gene_var)
    dim(sig_gene_var)
    table(sig_gene_var$pc)

    (num_genes <- length(unique(sig_gene_var$gene_id)))
    (num_features <- length(unique(sig_gene_var$feature_id)))
    (gwas_ld_blocks <- length(unique(sig_gene_var$gwas_ld_block)))
    (snp_summary <- summary(sig_gene_var$Num_SNP))

    final_out[['gene_3colocalized']] <- num_genes
    final_out[['feature_3colocalized']] <- num_features
    final_out[['gwas_ld_blocks_3colocalized']] <- gwas_ld_blocks
    final_out[['cluster_colocalized_snp_min']] <- snp_summary[[1]]
    final_out[['cluster_colocalized_snp_med']] <- snp_summary[[3]]
    final_out[['cluster_colocalized_snp_mean']] <- snp_summary[[4]]
    final_out[['cluster_colocalized_snp_max']] <- snp_summary[[6]]

    final_out %>%
        gather('count_type', 'count')
}

get_ptwas_counts <- function(runname, ptwas_thresh){
    # get ptwas counts
    datafile <- paste0("/path/to/new_rosmap/output/07_ptwas/07_ptwas_scan/wightman/with_proportions/",
                       runname, "/all_chroms_ptwas_scan.stratified_out.txt")
    ptwas <- read_table(datafile)
    head(ptwas)

    # filter for significant zscores
    sig_ptwas <- ptwas %>%
        filter(abs(stat) >= ptwas_thresh) %>%
        rename(feature_id = gene) %>%
        separate(feature_id, c('gene_id', 'pc'), sep="\\.PC", remove=FALSE, fill='right') 
    head(sig_ptwas)
    dim(sig_ptwas)

    (num_genes <- length(unique(sig_ptwas$gene_id)))
    (num_features <- length(unique(sig_ptwas$feature_id)))

    (final_out <- data.frame('gene_4ptwas'=num_genes,
                             'feature_4ptwas'=num_features))


    # check more detailed results
    datafile <- paste0("/path/to/new_rosmap/output/07_ptwas/07_ptwas_scan/wightman/with_proportions/",
                       runname, "/all_chroms_ptwas_scan.summary_out.txt")
    ptwas_summary <- read_table(datafile)
    head(ptwas_summary)
    dim(ptwas_summary)

    # filter for those found with high z score
    sig_summary <- ptwas_summary %>%
        filter(gene %in% sig_ptwas$feature_id)
    head(sig_summary)
    dim(sig_summary)

    # check out number of snps per gene
    (snp_summary <- summary(sig_summary$n_snps))

    final_out[['cluster_ptwas_min']] <- snp_summary[[1]]
    final_out[['cluster_ptwas_med']] <- snp_summary[[3]]
    final_out[['cluster_ptwas_mean']] <- snp_summary[[4]]
    final_out[['cluster_ptwas_max']] <- snp_summary[[6]]

    final_out %>%
        gather('count_type', 'count')
}

get_intact_counts <- function(runname, intact_thresh){
    # get intact counts
    datafile <- paste0("/path/to/new_rosmap/output/07_ptwas/11_intact/wightman/", 
                       runname, "_all_intact_results.tsv") 
    intact <- readRDS(datafile)
    head(intact)

    sig_intact <- intact %>%
        filter(intact_pip > intact_thresh) %>%
        separate(gene, c('gene_id', 'pc'), sep="\\.PC", fill='right')
    head(sig_intact)

    (num_genes <- length(unique(sig_intact$gene_id)))
    (num_features <- nrow(sig_intact))

    data.frame(count_type=c('gene_5intact', 'feature_5intact'), 
              count=c(num_genes, num_features) 
    )
}

plot_counts <- function(counts_df){
    head(counts_df)

    # format the counts table nicer
    formatted_df <- counts_df %>%
        gather('count_type', 'count', -summary_type, -cell_type, -region_type) 
        #separate(runname, c('summary_type', 'cell_type_region_type'), sep='_', extra='merge') %>%
        #mutate(cell_type_region_type = str_replace(cell_type_region_type, 'oligo_opc', 'OligoOPC')) %>%
        #separate(cell_type_region_type, c('cell_type', 'region_type'), sep='_', extra='merge')
    head(formatted_df)
    unique(formatted_df$summary_type)
    unique(formatted_df$cell_type)
    unique(formatted_df$region_type)

    # save a wide version of the plot
    formatted_wide <- formatted_df %>%
        spread(count_type, count)
    head(formatted_wide)

    (savefile <- paste0(savedir, "formatted_counts_wide.csv"))
    write_csv(formatted_wide, savefile)

    # plot proportions
    head(formatted_df)
    keep_rt <- 'full_gene'
    keep_counts <- c('gene_0tested', 'gene_1qtl', 'gene_2finemapped')
    p <- formatted_df %>%
        mutate(count = as.numeric(count)) %>%
        left_join(cell_type_map) %>%
        left_join(summary_type_map) %>%
        filter(
               #count < 1e5,
               region_type == keep_rt,
               count_type %in% keep_counts) %>%
        #spread(count_type, count) %>%
        #mutate(qtl_prop = round(gene_1qtl / gene_0tested * 100, 2)) %>%
        #mutate(finemapped_prop = round(gene_2finemapped / gene_0tested * 100, 2)) %>%
        #select(summary_type:formatted_st, qtl_prop, finemapped_prop) %>%
        #gather('prop_type', 'proportion', qtl_prop:finemapped_prop) %>%
        #mutate(prop_type = fct_recode(prop_type,
                                       #'QTL mapping'='qtl_prop',
                                       #'Fine-mapping'='finemapped_prop'
               #)) %>%
        #mutate(prop_type = fct_relevel(prop_type, 'QTL mapping', 'Fine-mapping')) %>%
        filter(count_type %in% c('gene_2finemapped', 'gene_1qtl')) %>%
        mutate(prop_type= fct_recode(count_type,
                                       'Fine-mapping'='gene_2finemapped',
                                       'QTL mapping'='gene_1qtl'
               )) %>%
        mutate(prop_type = fct_relevel(prop_type, 'QTL mapping', 'Fine-mapping')) %>%
        ggplot(aes(x=formatted_st, y=count, 
                   #group=interaction(formatted_ct, formatted_st), 
                   color=formatted_ct
                   #shape=formatted_st
               )
        ) +
        geom_point(size=3) +
        geom_line(aes(group=formatted_ct), linetype='dashed') +
        theme_bw() +
        facet_wrap(. ~ prop_type, scales='free') +
        theme_bw() +
        theme(text = element_text(size=20),
              axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5),
              strip.background = element_rect(fill='white', colour='white'),
              strip.text = element_text(color='black', face='bold'),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
              ) +
        scale_color_manual(values=cell_type_colors, name="Cell type") +
        #labs(linetype="Summary type") +
        #labs(shape="Summary type") +
        #guides(linetype=guide_legend(keywidth = 2.5, keyheight = 1, name="Summary type")) +
        ylab(paste("Significant meGene count")) +
        xlab(paste("Summary type")) 
    #(savefile <- paste0(plotdir, keep_rt, "_proportion_plot.png"))
    (savefile <- paste0(plotdir, keep_rt, "_count_simple_plot.png"))
    ggsave(p, file=savefile, height=5, width=8)
    

    # plot the counts
    head(formatted_df)
    unique(formatted_df$count_type)
    keep_counts <- c('gene_0tested', 'gene_1qtl', 'gene_2finemapped')
    keep_rt <- 'promoters'
    p1 <- formatted_df %>%
        mutate(count = as.numeric(count)) %>%
        left_join(cell_type_map) %>%
        left_join(summary_type_map) %>%
        filter(
               #count < 1e5,
               region_type == keep_rt,
               count_type %in% keep_counts) %>%
        mutate(count_type = fct_recode(count_type,
                                       'Fine-mapping'='gene_2finemapped',
                                       'QTL mapping'='gene_1qtl',
                                       'Total'='gene_0tested'
               )) %>%
        ggplot(aes(x=summary_type, y=count, 
                   group=interaction(formatted_ct, formatted_st), 
                   color=formatted_ct, linetype=formatted_st, shape=formatted_st)) +
        #geom_line(linewidth=2) +
        #geom_point(size=4, alpha=1) +
        geom_jitter(size=4, alpha=1, height=0) +
        #scale_y_continuous(trans='log10') +
        facet_wrap(. ~ count_type, scales='free') +
        theme_bw() +
        theme(text = element_text(size=20),
              axis.text.x = element_text(angle=20, hjust=1, vjust=1),
              strip.background = element_rect(fill='white', colour='white'),
              strip.text = element_text(color='black', face='bold'),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
              ) +
        scale_color_manual(values=cell_type_colors, name="Cell type") +
        #labs(linetype="Summary type") +
        labs(shape="Summary type") +
        guides(linetype=guide_legend(keywidth = 2.5, keyheight = 1, name="Summary type")) +
        ylab(paste("Gene count")) +
        xlab(paste("Steps in QTL analysis")) 


    p2 <- formatted_df %>%
        mutate(count = as.numeric(count)) %>%
        left_join(cell_type_map) %>%
        left_join(summary_type_map) %>%
        filter(
               #count < 1e5,
               region_type == keep_rt,
               count_type %in% keep_counts) %>%
        mutate(count_type = fct_recode(count_type,
                                       'Fine-mapped mGenes'='gene_2finemapped',
                                       'QTL mGenes'='gene_1qtl',
                                       'Total'='gene_0tested'
               )) %>%
        ggplot(aes(x=count_type, y=count, 
                   group=interaction(formatted_ct, formatted_st), 
                   color=formatted_ct, linetype=formatted_st, shape=formatted_st)) +
        geom_line(linewidth=2) +
        geom_point(size=4, alpha=1) +
        #geom_jitter(size=4, alpha=1, height=0) +
        #scale_y_continuous(trans='log10') +
        #facet_wrap(. ~ count_type, scales='free') +
        theme_bw() +
        theme(text = element_text(size=20),
              axis.text.x = element_blank(),
              #axis.text.x = element_text(angle=20, hjust=1, vjust=1),
              axis.title.x = element_blank(),
              strip.background = element_rect(fill='white', colour='white'),
              strip.text = element_text(color='black', face='bold'),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.box.margin = margin(r=20)
              ) +
        scale_color_manual(values=cell_type_colors, name="Cell type") +
        #labs(linetype="Summary type") +
        labs(shape="Summary type", linetype="Summary type") +
        guides(linetype=guide_legend(keywidth = 2.5, keyheight = 1, name="Summary type")) +
        ylab(paste("Gene count"))
        #xlab(paste("Steps in QTL analysis")) 


    combined_p <- ggarrange(p2, p1, 
                            nrow=2, ncol=1,
                            common.legend =TRUE,
                            legend='right'
    )

    (savefile <- paste0(plotdir, keep_rt, "_counts_line_point.png"))
    ggsave(combined_p, file=savefile, width=10, height=11)
}

main <- function(){
    summary_types <- c('avgs', 'pcs')
    region_types <- c('full_gene', 'promoters', 'preTSS', 'gene_body')
    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')

    # get all combinations of parameters
    runs <- expand.grid(summary_type = summary_types,
                        region_type = region_types,
                        cell_type = cell_types
    ) %>%
        mutate(rn = row_number())
    head(runs)

    # define thresholds for significance
    coloc_thresh = 0.5
    intact_thresh = 0.5
    ptwas_thresh = 5.45 # genome wide significance


    # process each run
    row <- runs[1,]
    process_run <- function(row){
        summary_type <- row[['summary_type']]
        region_type <- row[['region_type']]
        cell_type <- row[['cell_type']]
        rn <- row[['rn']]

        (runname <- paste(summary_type, cell_type, region_type, sep='_'))
        print(paste("Processing counts for", runname, rn, "out of", nrow(runs), Sys.time()))

        # get counts of total gene x variants tested
        print("Getting total counts")
        total_counts <- get_total_counts(runname)
        total_counts

        # get counts of sig QTLs from QTL mapping
        print("Getting significant QTL counts")
        old_runname <- paste(cell_type, region_type, summary_type, sep='_')
        sig_counts <- get_sig_mapped_qtls(old_runname)

        # get counts of sig fine-mapped QTLs
        print("Getting finemapped counts")
        finemapped_counts <- get_finemapped_qtls(runname)

        # get counts of colocalizations
        print("Getting colocalization counts")
        (colocalization_counts <- get_colocalization_counts(runname, coloc_thresh))

        # get counts of PTWAS
        print("Getting PTWAS counts")
        (ptwas_counts <- get_ptwas_counts(runname, ptwas_thresh))

        # get counts of INTACT
        print("Getting INTACT counts")
        (intact_counts <- get_intact_counts(runname, intact_thresh))

        all_counts <- rbind(total_counts, finemapped_counts) %>%
            rbind(., sig_counts) %>%
            rbind(., colocalization_counts) %>%
            rbind(., ptwas_counts) %>%
            rbind(., intact_counts) %>%
            spread(count_type, count) %>%
            mutate(summary_type = summary_type,
                   cell_type = cell_type,
                   region_type = region_type
            )

        (savefile <- paste0(savedir, runname, "_result_counts.tsv"))
        write_tsv(all_counts, savefile)

        all_counts
    }

    counts <- apply(runs, 1, process_run)
    counts_df <- do.call(rbind, counts)
    head(counts_df)
    dim(counts_df)

    savefile <- paste0(savedir, "all_counts.tsv")
    write_tsv(counts_df, savefile)

    plotting <- FALSE
    if (plotting){
        counts_df <- read_tsv(savefile)
        plot_counts(counts_df)
    }
}

main()

