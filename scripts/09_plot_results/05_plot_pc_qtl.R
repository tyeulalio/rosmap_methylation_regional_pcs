library(ggplot2) 
library(genio)
library(tidyverse)

# plot qtl at the variants x methylation -level to visualize
# a qtl that was identified as sig by PCs but not averages

savedir <- paste0("/path/to/output/")
dir.create(savedir)

load_qtl_list <- function(){
    # load the pre-filtered list of 
    # qtls of interest (PCs but not averages)
    datafile <- paste0("/path/to/new_rosmap/output/05_qtl_analysis/13_check_qtl_results/with_proportions/sigPC_nosigAvgs_qtls.csv")
    qtl_list <- read_csv(datafile)
    head(qtl_list)
    qtl_list
}

select_qtl <- function(qtl_list){
    # select the qtl to plot
    # go simple for now and choose top one
    head(qtl_list)

    (top_genes <- unique(qtl_list$symbol) %>% head())
    filtered_list <- qtl_list %>%
        filter(symbol %in% top_genes)

    (row <- sample(1:nrow(filtered_list), 1))

    qtl <- filtered_list[row,]
    qtl
}

load_variant_info <- function(qtl){
    # load the variant info
    qtl

    # grab variant ID
    variant_id <- qtl[['variant_id']]
    variant_id

    # read annotation tables
    data_prefix <- paste0("/path/to/new_rosmap/output/05_qtl_analysis/09_remove_outliers/all_chroms_hg38")
    datafile <- paste0(data_prefix, ".bim")
    bim <- read_bim(datafile)

    # read fam
    datafile <- paste0(data_prefix, ".fam")
    fam <- read_fam(datafile)
    
    # read bed
    datafile <- paste0(data_prefix, ".bed")
    bed <- read_bed(datafile, bim$id, fam$id)
    head(bed)

    # find variant
    variant_info <- bed[variant_id,] %>%
        as.data.frame() %>%
        rename(genotype=".") %>%
        rownames_to_column('subject_id')
    head(variant_info)

    variant_info
}

load_methylation_info <- function(qtl){
    # load methylation data 
    qtl

    # grab relevant info
    cell_type <- qtl[['cell_type']]
    region_type <- qtl[['region_type']]
    qtl_gene <- qtl[['gene_id']]

    summary_type <- 'pcs'
    load_summary_type_meth <- function(summary_type){
        filename <- paste(cell_type, region_type, "formatted", summary_type, sep='_')
        datafile <- paste0("/path/to/new_rosmap/output/05_qtl_analysis/03_format_pc_beds/",
                           filename, ".bed.gz")
        st_meth <- read_table(datafile)
        head(st_meth)

        gene_meth <- st_meth %>%
            filter(grepl(qtl_gene, gene_id)) %>%
            select(-"#chr", -start, -end) %>%
            gather('subject_id', 'meth', -gene_id) %>%
            mutate(summary_type = summary_type) %>%
            separate(gene_id, c('gene_id', 'pc'), sep='-', fill='right')
        gene_meth
    }
    summary_types <- c('avgs', 'pcs')
    qtl_meth <- lapply(summary_types, load_summary_type_meth)
    qtl_meth_df <- do.call(rbind, qtl_meth)
    head(qtl_meth_df)

    qtl_meth_df
}


plot_info <- function(variant_info, meth_info, qtl){
    # plot the variant and methylation info together
    head(variant_info)
    head(meth_info)

    # join the data togethr
    joined_info <- variant_info %>%
        left_join(meth_info)
    head(joined_info)

    # check if there are multiple PCs
    unique(joined_info$pc) 
    if (length(unique(joined_info$pc)) > 2) stop("Multiple PCs for this gene")

    gene <- qtl[['symbol']]
    variant <- qtl[['variant_id']]
    variant

    formatted_variant <- variant %>%
        str_split('_') %>%
        .[[1]] %>%
        .[1:2] %>%
        str_c(collapse=":")
    formatted_variant

    # plot the data
    p <- ggplot(joined_info) +
        geom_boxplot(aes(x=as.factor(genotype), y=meth, fill=summary_type)) +
        facet_wrap(. ~ summary_type, scales='free') +
        theme_bw() +
        ylab(paste(gene, "normalized methylation")) +
        xlab(paste("Genotype (", formatted_variant, ")"))
    (savefile <- paste0(savedir, gene, "_", variant, "_qtl_variant_meth_boxplot.png"))
    ggsave(p, file=savefile, height=5, width=6)

    # save the plot object
    savefile <- paste0(savedir, gene, "_", variant, "_qtl_variant_meth_plot.rds")
    saveRDS(p, savefile)
    1
}

main <- function(){
    # load the qtl list
    qtl_list <- load_qtl_list()
    head(qtl_list)

    set.seed(1174117174)
    for (i in 1:10){
        # select the qtl to focus on
        qtl <- select_qtl(qtl_list)
        print(qtl)

        # load variant info 
        variant_info <- load_variant_info(qtl)

        # load methylation info
        meth_info <- load_methylation_info(qtl)

        # plot the variant and methylation info together
        plot_info(variant_info, meth_info, qtl)
    }

}
