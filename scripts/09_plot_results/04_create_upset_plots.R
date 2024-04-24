# create upset plots
library(ggplot2)
library(ggplotify)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)
#library(cowplot)
library(ComplexUpset)
library(UpSetR)
library(tidyverse)

savedir <- paste0("/path/to/output/")
dir.create(savedir)


get_total_features <- function(runname){
    # get total genes tested
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/01_torus_input/with_proportions/", runname, "_fastqtl_single_snp_output.tsv")

    # get all features
    (cmd <- paste0("cut -d' ' -f1 ", datafile, " | uniq"))
    features <- system(cmd, intern=TRUE)
    head(features)
    length(features)

    features
}

# load the qtl genes 
get_qtl_genes <- function(cell_type, region_type, summary_type){
    old_runname <- paste(cell_type, region_type, summary_type, sep='_')
    print(paste("Loading qtl data for", old_runname))

    # get sig qtls after qtl mapping
    datafile <- paste0("/path/to/new_rosmap/output/05_qtl_analysis/12_map_qtls/with_proportions/",
                       old_runname, "/cis_qtl.signif_pairs.csv")
    sig_qtls <- read_csv(datafile) %>%
        separate(phenotype_id, c('phenotype_id', 'pc'), sep='-', fill='right')
    head(sig_qtls)

    # get sig qtl genes
    qtl_genes <- sig_qtls %>%
        select(phenotype_id) %>%
        unique() %>%
        separate(phenotype_id, c('gene', 'gene_vers'), sep='\\.') 
    head(qtl_genes)

    unique(qtl_genes$gene)
}

get_finemapped_genes <- function(cell_type, region_type, summary_type){
    runname <- paste(summary_type, cell_type, region_type, sep='_')
    print(paste("Loading qtl data for", runname))

    # get counts of significant finemapped qtls
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/05_fastqtl_qtl_annotations/with_proportions/", runname, "/fastenloc.qtl.annotation.vcf.gz")
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


    # count unique genes
    finemapped_genes <- formatted_finemapped %>%
        select(gene_id) %>%
        unique() %>%
        separate(gene_id, c('gene', 'gene_vers'), sep='-')
    head(finemapped_genes)

    unique(finemapped_genes$gene)
}

get_coloc_genes <- function(cell_type, region_type, summary_type){
    runname <- paste(summary_type, cell_type, region_type, sep='_')
    print(paste("Loading qtl data for", runname))

    # get counts of colocalizations
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/04_fastenloc/01_fastenloc_output/with_proportions/", runname, "/wightman/fastqtl_ld.enloc.gene.out")
    gene_colocalizations <- read_table(datafile)
    head(gene_colocalizations)

    # filter for high colocalization probabilities
    sig_colocalizations <- gene_colocalizations %>%
        filter(GLCP > coloc_thresh) %>%
        separate(Gene, c('gene_id', 'pc'), sep='-PC', fill='right') %>%
        separate(gene_id, c('gene', 'gene_vers'), sep='-')
    sig_colocalizations

    unique(sig_colocalizations$gene)
}

get_ptwas_genes <- function(cell_type, region_type, summary_type){
    runname <- paste(summary_type, cell_type, region_type, sep='_')
    print(paste("Loading qtl data for", runname))

    # get ptwas counts
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/07_ptwas_scan/wightman/with_proportions/",
                       runname, "/all_chroms_ptwas_scan.stratified_out.txt")
    ptwas <- read_table(datafile)
    head(ptwas)

    # filter for significant zscores
    sig_ptwas <- ptwas %>%
        filter(abs(stat) >= ptwas_thresh) %>%
        rename(feature_id = gene) %>%
        separate(feature_id, c('gene_id', 'pc'), sep="\\.PC", remove=FALSE, fill='right') %>%
        select(gene_id) %>%
        separate(gene_id, c('gene', 'gene_vers'), sep='\\.')
    head(sig_ptwas)
    dim(sig_ptwas)

    unique(sig_ptwas$gene)
}

get_intact_genes <- function(cell_type, region_type, summary_type){
    runname <- paste(summary_type, cell_type, region_type, sep='_')
    print(paste("Loading qtl data for", runname))

    # get intact counts
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/11_intact/wightman/", 
                       runname, "_all_intact_results.tsv") 
    intact <- readRDS(datafile)
    head(intact)

    sig_intact <- intact %>%
        filter(intact_pip > intact_thresh) %>%
        separate(gene, c('gene_id', 'pc'), sep="\\.PC", fill='right') %>%
        select(gene_id) %>%
        separate(gene_id, c('gene', 'gene_novers'), sep='\\.')
    head(sig_intact)

    unique(sig_intact$gene)
}

# checkign the qtl resuls
check_qtls <- function(cell_types, region_type, summary_type, colors_map){ 
    # load genes for each cell type
    cell_type <- 'astro'
    load_data <- function(cell_type){
        get_qtl_genes(cell_type, region_type, summary_type)
    }
    qtl_genes <- lapply(colors_map$cell_type, load_data) 

    # create cell type color map
    ct_map <- data.frame(cell_type=cell_types) %>%
        left_join(cell_type_map)
    names(qtl_genes) <- ct_map$formatted_ct

    # create an upset plot using upsetr

    # create an upset plot using upsetr
    (savefile <- paste0(savedir, summary_type, "_", region_type,
                       "_celltype_0qtl_genes_overlaps_upset"
    ))
    angle=0
    width=10
    height=8
    p <- plot_upset(qtl_genes, savefile, colors_map, angle, width, height)

    p
}

#gene_lists <- colocalized_genes
#genes_list <-finemapped_genes
plot_upset <- function(gene_lists, savefile, colors_map, angle, width, height){
    # need to get order of gene list lengths
    gene_lengths <- data.frame(cell_type=colors_map$cell_type, set_length=unlist(lapply(gene_lists, length))) %>%
        arrange(-set_length) %>%
        left_join(colors_map) 
    gene_lengths

    # set the colors correctly
    sets_color <- gene_lengths$cell_type_colors
    names(sets_color) <- gene_lengths$formatted_ct
    sets_color


    # plot upset plot
    pdf_savefile <- paste0(savefile, ".pdf")
    pdf(pdf_savefile, width=width, height=8)
    p <- upset(fromList(gene_lists),
               number.angles=angle,
               order.by='freq',
               sets.bar.color=sets_color,
               nintersects=10,
               text.scale=c(3,3,3,2,3,3)
    )
    print(p)
    dev.off()

    rds_savefile <- paste0(savefile, ".rds")
    saveRDS(p, file=rds_savefile)


    return(p)
}

check_finemapped <- function(cell_types, region_type, summary_type, colors_map){ 
    # load genes for each cell type
    cell_type <- 'astro'
    load_data <- function(cell_type){
        runname <- paste(summary_type, cell_type, region_type, sep='_')
        print(paste("Loading qtl data for", runname))

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


        # count unique genes
        finemapped_genes <- formatted_finemapped %>%
            select(gene_id) %>%
            unique() %>%
            separate(gene_id, c('gene', 'gene_vers'), sep='-')
        head(finemapped_genes)

        unique(finemapped_genes$gene)
    }
    finemapped_genes <- lapply(cell_types, load_data) 
    ct_map <- data.frame(cell_type=cell_types) %>%
        left_join(cell_type_map)
    head(ct_map)

    names(finemapped_genes) <- ct_map$formatted_ct

    # create an upset plot using upsetr
    (savefile <- paste0(savedir, summary_type, "_", region_type,
                       "_celltype_1finemapped_genes_overlaps_upset"
    ))
    angle=0
    width=10
    height=8
    p <- plot_upset(finemapped_genes, savefile, colors_map, angle, width, height)

    return(p)
}

check_colocalization <- function(cell_types, region_type, summary_type, coloc_thresh, colors_map){ 
    # load genes for each cell type
    cell_type <- 'astro'
    load_data <- function(cell_type){
        runname <- paste(summary_type, cell_type, region_type, sep='_')
        print(paste("Loading qtl data for", runname))

        # get counts of colocalizations
        datafile <- paste0("/path/to/new_rosmap/output/06_fine_mapping_colocalization/04_fastenloc/01_fastenloc_output/with_proportions/", runname, "/wightman/fastqtl_ld.enloc.gene.out")
        gene_colocalizations <- read_table(datafile)
        head(gene_colocalizations)

        # filter for high colocalization probabilities
        sig_colocalizations <- gene_colocalizations %>%
            filter(GLCP > coloc_thresh) %>%
            separate(Gene, c('gene_id', 'pc'), sep='-PC', fill='right') %>%
            separate(gene_id, c('gene', 'gene_vers'), sep='-')
        sig_colocalizations

        unique(sig_colocalizations$gene)
    }
    colocalized_genes <- lapply(cell_types, load_data) 
    ct_map <- data.frame(cell_type=cell_types) %>%
        left_join(cell_type_map)
    head(ct_map)

    names(colocalized_genes) <- ct_map$formatted_ct

    head(colocalized_genes)

    # get the overlaps between groups
    gene_list <- colocalized_genes

    # create an upset plot using upsetr
    savefile <- paste0(savedir, summary_type, "_", region_type,
                       "_celltype_2colocalized_genes_overlaps_upset"
    )
    savefile
    p <- plot_upset(colocalized_genes, savefile, colors_map)

    return(p)
}

get_overlap_stats <- function(genes_list){
    all_genes <- stack(genes_list)
    colnames(all_genes) <- c('gene_id', 'cell_type')
    head(all_genes)

    # create indicator dataframe for each gene across cell types
    genes_df <- all_genes %>%
        mutate(exists = 1) %>%
        spread(cell_type, exists, fill=0)
    head(genes_df)

    # store cell types
    cts <- colnames(genes_df)[2:ncol(genes_df)]

    # convert group indicator columns to boolean
    genes_df[cts] <- genes_df[cts] == 1

    # get intersection sizes
    ints <- data.frame(expected_count=upset_data(genes_df, intersect=cts)$sizes$exclusive_intersection) %>%
        arrange(-expected_count) %>%
        rownames_to_column('intersect') %>%
        mutate(cell_type_specific = intersect %in% cts) 
    head(ints)
    ints

    summarized_ints <- ints %>%
        group_by(cell_type_specific) %>%
        summarize(count = sum(expected_count))
    head(summarized_ints)

    print(ints)
}

check_ptwas <- function(cell_types, region_type, summary_type, ptwas_thresh, colors_map){ 
    # load genes for each cell type
    cell_type <- 'astro'
    load_data <- function(cell_type){
        runname <- paste(summary_type, cell_type, region_type, sep='_')
        print(paste("Loading qtl data for", runname))

        # get ptwas counts
        datafile <- paste0("/path/to/new_rosmap/output/07_ptwas/07_ptwas_scan/wightman/with_proportions/",
                           runname, "/all_chroms_ptwas_scan.stratified_out.txt")
        ptwas <- read_table(datafile)
        head(ptwas)

        # filter for significant zscores
        sig_ptwas <- ptwas %>%
            filter(abs(stat) >= ptwas_thresh) %>%
            rename(feature_id = gene) %>%
            separate(feature_id, c('gene_id', 'pc'), sep="\\.PC", remove=FALSE, fill='right') %>%
            select(gene_id) %>%
            separate(gene_id, c('gene', 'gene_vers'), sep='\\.')
        head(sig_ptwas)
        dim(sig_ptwas)

        unique(sig_ptwas$gene)
    }
    ptwas_genes <- lapply(cell_types, load_data) 
    names(ptwas_genes) <- cell_types

    get_overlap_stats(ptwas_genes)

    # create an upset plot using upsetr
    savefile <- paste0(savedir, summary_type, "_", region_type,
                       "_celltype_3ptwas_genes_overlaps_upset"
    )
    p <- plot_upset(ptwas_genes, savefile, colors_map)

    return(p)
}

check_intact <- function(cell_types, region_type, summary_type, intact_thresh, colors_map){ 
    # load genes for each cell type
    cell_type <- 'astro'
    load_data <- function(cell_type){
        runname <- paste(summary_type, cell_type, region_type, sep='_')
        print(paste("Loading qtl data for", runname))

        # get intact counts
        datafile <- paste0("/path/to/new_rosmap/output/07_ptwas/11_intact/wightman/", 
                           runname, "_all_intact_results.tsv") 
        intact <- readRDS(datafile)
        head(intact)

        sig_intact <- intact %>%
            filter(intact_pip > intact_thresh) %>%
            separate(gene, c('gene_id', 'pc'), sep="\\.PC", fill='right') %>%
            select(gene_id) %>%
            separate(gene_id, c('gene', 'gene_novers'), sep='\\.')
        head(sig_intact)

        unique(sig_intact$gene)
    }
    intact_genes <- lapply(cell_types, load_data) 
    names(intact_genes) <- cell_types

    # create an upset plot using upsetr
    savefile <- paste0(savedir, summary_type, "_", region_type,
                       "_celltype_4intact_genes_overlaps_upset"
    )
    p <- plot_upset(intact_genes, savefile, colors_map)

    return(p)
}

combine_plots <- function(p0, p1, p2, p3, p4, region_type, summary_type, colors_map){
    # combine the plots into one
    p <- p2
    format_plot <- function(p, title){
        formatted_p <- as.ggplot(p) +
            ggtitle(title)
        formatted_p
    }

    formatted_p0 <- format_plot(p0, 'QTL mapping')
    formatted_p1 <- format_plot(p1, 'Fine-mapping')
    formatted_p2 <- format_plot(p2, 'Colocalization')
    formatted_p3 <- format_plot(p3, 'PTWAS')
    formatted_p4 <- format_plot(p4, 'INTACT')

    combined_plots <- grid.arrange(
                                   formatted_p0,
                                   formatted_p1,
                                   formatted_p2,
                                   formatted_p3,
                                   formatted_p4,
                                   nrow=2,
                                   top=paste(summary_type, region_type, "Gene Counts")
    )
    plotdir <- paste0(savedir, "combined_plots/")
    dir.create(plotdir, showWarnings = FALSE)
    (savefile <- paste0(plotdir, summary_type, "_", region_type, "_combined_genes_overlap_upset.png"))
    ggsave(combined_plots, file=savefile, height=10, width=16)
    1
}

gene_lists <- genes
plot_upset_general <- function(gene_lists, savefile){
    # need to get order of gene list lengths
    #gene_lengths <- data.frame(cell_type=colors_map$cell_type, set_length=unlist(lapply(gene_lists, length))) %>%
        #arrange(-set_length) %>%
        #left_join(colors_map) 
    #gene_lengths

    # set the colors correctly
    names(gene_lists)
    summary_type_colors
    sets_color <- summary_type_colors[names(gene_lists)]
    #names(sets_color) <- gene_lengths$formatted_ct
    sets_color


    # plot upset plot
    pdf(savefile, width=6, height=6)
    p <- upset(fromList(gene_lists),
               order.by='freq',
               sets.bar.color=sets_color,
               nintersects=10,
               text.scale=c(2,2,2,1.25,2,2)
               #mb.ratio = c(0.8,0.2)
    )
    print(p)
    dev.off()

    return(p)
}

plots_list <- ct_plots
combine_plots_st <- function(plots_list, region_type, qtl_step){
    # combine the plots into one
    format_plot <- function(p, title){
        formatted_p <- as.ggplot(p) +
            ggtitle(title)
        formatted_p
    }

    names(plots_list)
    names(plots_list) <- c("Astrocyte", "Endothelial", "Neuron", "Oligo/OPC", "Bulk")

    i <- 1
    formatted_p1 <- format_plot(plots_list[[i]], names(plots_list)[[i]])
    i <- 2
    formatted_p2 <- format_plot(plots_list[[i]], names(plots_list)[[i]])
    i <- 3
    formatted_p3 <- format_plot(plots_list[[i]], names(plots_list)[[i]])
    i <- 4
    formatted_p4 <- format_plot(plots_list[[i]], names(plots_list)[[i]])
    i <- 5
    formatted_p5 <- format_plot(plots_list[[i]], names(plots_list)[[i]])

    combined_plots <- grid.arrange(formatted_p1,
                                formatted_p2,
                                formatted_p3,
                                formatted_p4,
                                formatted_p5,
                                   nrow=3,
                                   top=paste(qtl_step, region_type, "Gene Counts")
    )

    plotdir <- paste0(savedir, "combined_plots/")
    dir.create(plotdir, showWarnings = FALSE)
    (savefile <- paste0(plotdir, "summary_type_", qtl_step, "_", region_type, "_combined_genes_overlap_upset.png"))
    ggsave(combined_plots, file=savefile, height=10, width=8)
    1
}

region_type <- 'full_gene'
compare_summary_types <- function(region_type){
    # make umaps to compare summary types

    cell_type <- cell_types[[1]]
    get_genes <- get_finemapped_genes
    qtl_step <- 'finemapping'
    get_ct_plots <- function(cell_type, region_type, get_genes, qtl_step){
        print(paste("Loading data for", cell_type, region_type))
        # load the qtl data
        load_data <- function(summary_type){ 
            genes <- get_genes(cell_type, region_type, summary_type)
            genes
        }
        genes <- lapply(summary_types, load_data)
        names(genes) <- c('Avgs', 'rPCs')
        head(genes)
        
        (savefile <- paste0(savedir, cell_type, "_", region_type, 
                           "_summary_type_", qtl_step, "_upset.pdf"))
        p <- plot_upset_general(genes, savefile)
        p
    }

    print("Getting QTL counts")
    # put all cell types in one plot
    #get_ct_qtl_plots <- function(cell_type, region_type){
        #print(paste("Loading data for", cell_type, region_type))
        ## load the qtl data
        #load_data <- function(summary_type){ 
            #genes <- get_qtl_genes(cell_type, region_type, summary_type)
            #genes
        #}
        #genes <- lapply(summary_types, load_data)
        #names(genes) <- summary_types
        
        #p <- plot_upset_general(genes)
        #p
    #}
    qtl_step <- "qtl_mapping"
    ct_plots <- lapply(cell_types, get_ct_plots, region_type, get_qtl_genes, qtl_step)
    names(ct_plots) <- cell_types
    combined_p <- combine_plots_st(ct_plots, region_type, qtl_step)


    # get counts of sig fine-mapped QTLs
    print("Getting finemapped counts")
    qtl_step <- "finemapping"
    ct_plots <- lapply(cell_types, get_ct_plots, region_type, get_finemapped_genes, qtl_step)
    names(ct_plots) <- cell_types
    combined_p <- combine_plots_st(ct_plots, region_type, qtl_step)


    # get counts of colocalizations
    print("Getting colocalization counts")
    ct_plots <- lapply(cell_types, get_ct_plots, region_type, get_coloc_genes)
    names(ct_plots) <- cell_types
    qtl_step <- "colocalization"
    combined_p <- combine_plots_st(ct_plots, region_type, qtl_step)

    # get counts of PTWAS
    print("Getting PTWAS counts")
    ct_plots <- lapply(cell_types, get_ct_plots, region_type, get_ptwas_genes)
    names(ct_plots) <- cell_types
    qtl_step <- "ptwas"
    combined_p <- combine_plots_st(ct_plots, region_type, qtl_step)

    # get counts of INTACT
    print("Getting INTACT counts")
    ct_plots <- lapply(cell_types, get_ct_plots, region_type, get_intact_genes)
    names(ct_plots) <- cell_types
    qtl_step <- "ptwas"
    combined_p <- combine_plots_st(ct_plots, region_type, qtl_step)

    1
}

load_celltype_plot <- function(region_type, summary_type, method){
    savefile <- paste0(savedir, summary_type, "_", region_type,
                       "_celltype_", method, "_genes_overlaps_upset.rds"
    )
    p <- readRDS(savefile)
    p
}

combine_summary_type_plots <- function(runs){
    # make qtl-mapping/fine-mapping plot
    region_type <- 'full_gene'
    combined_qtl_mapping <- function(region_type){
        # load qtl plots
        method <- '0qtl'
        summary_type <- 'pcs'
        pc_qtl <- load_celltype_plot(region_type, summary_type, method)
        summary_type <- 'avgs'
        avg_qtl <- load_celltype_plot(region_type, summary_type, method)

        # load fine-mapping plots
        method <- '1finemapped'
        summary_type <- 'pcs'
        pc_finemapped <- load_celltype_plot(region_type, summary_type, method)
        summary_type <- 'avgs'
        avg_finemapped <- load_celltype_plot(region_type, summary_type, method)

        # combine these plots together
        p <- pc_qtl
        format_plot <- function(p, title){
            formatted_p <- as.ggplot(p) +
                ggtitle(title) 
            formatted_p
        }

        formatted_p0 <- format_plot(avg_qtl, 'QTL mapping with averages')
        formatted_p1 <- format_plot(pc_qtl, 'QTL mapping with rPCs')
        formatted_p2 <- format_plot(avg_finemapped, 'Fine-mapping with averages')
        formatted_p3 <- format_plot(pc_finemapped, 'Fine-mapping with rPCs')

        combined_plots <- grid.arrange(
                                       formatted_p0,
                                       formatted_p1,
                                       formatted_p2,
                                       formatted_p3,
                                       nrow=2,
                                       top=paste(region_type, "QTL mapping")
        )
        plotdir <- paste0(savedir, "combined_plots/")
        dir.create(plotdir, showWarnings = FALSE)
        (savefile <- paste0(plotdir, region_type, "_combined_summarytypes_qtl_genes_overlap_upset.png"))
        ggsave(combined_plots, file=savefile, height=10, width=13)
        1

    }
}

main <- function(){
    summary_types <- c('avgs', 'pcs')
    region_types <- c('full_gene', 'promoters', 'preTSS', 'gene_body')
    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')


    coloc_thresh = 0.5
    intact_thresh = 0.5
    ptwas_thresh = 5.45 # genome wide significance

    # for testing
    summary_type <- 'avgs'
    region_type <- 'full_gene'

    # compare cell types in upset plots
    compare_celltypes <- function(region_type, summary_type, cell_types){
        # create a color map - these come from another file, must be sourced, it's colors for each cell type/summary type/region type
        cell_type_colors
        cell_type_map
        colors_df <- data.frame(cell_type_colors) %>%
            rownames_to_column('formatted_ct')
        colors_map <- cell_type_map %>%
            left_join(colors_df)
        colors_map

        # check sig qtls
        print("Checking significant QTLs")
        p0 <- check_qtls(cell_types, region_type, summary_type, colors_map)

        # get counts of sig fine-mapped QTLs
        print("Getting finemapped counts")
        p1 <- check_finemapped(cell_types, region_type, summary_type, colors_map)

        # get counts of colocalizations
        print("Getting colocalization counts")
        p2 <- check_colocalization(cell_types, region_type, summary_type, coloc_thresh, colors_map)

        # get counts of PTWAS
        print("Getting PTWAS counts")
        p3 <- check_ptwas(cell_types, region_type, summary_type, ptwas_thresh, colors_map)

        # get counts of INTACT
        print("Getting INTACT counts")
        p4 <- check_intact(cell_types, region_type, summary_type, intact_thresh, colors_map)

        combine_plots(p0, p1, p2, p3, p4, region_type, summary_type, colors_map)
        1
    }

    lapply(region_types, compare_celltypes, 'pcs', cell_types)
    lapply(region_types, compare_celltypes, 'avgs', cell_types)

    # combine the summary types into a single plot
    combine_summary_type_plots()
    

    
    region_type <- 'full_gene'
    cell_type <- 'astro'

    
    compare_summary_types(region_type, cell_type)
}

main()

