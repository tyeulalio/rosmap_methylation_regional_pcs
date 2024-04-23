library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

# combine result from coloc, twas, intact

source("/home/eulalio/deconvolution/new_rosmap/scripts/shared_functions/02_plot_colors.R")

gwas_type <- 'wightman'

(savedir <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/12_combined_results/", gwas_type, "/"))
dir.create(savedir)

load_data <- function(runs){
    # load data from the different analyses
    row <- runs[1,]

    # load intact results
    load_intact <- function(row){
        ct <- row[['cell_type']] %>% as.character()
        rt <- row[['region_type']] %>% as.character()
        rn <- row[['rn']] %>% as.character()
        st <- row[['summary_type']] %>% as.character()
        runname <- row[['runname']] %>% as.character()

        print(paste("Loading data for", ct, rt, st, rn, "out of", nrow(runs)))

        #print(paste("LOADING OLD DATA WITHOUT GLOBAL PC CLEANING!!!!"))
        #datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/without_global_pc_cleaning/11_intact/", gwas_type, "/", runname, "_all_intact_results.tsv")
        datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/11_intact/", gwas_type, "/", runname, "_all_intact_results.tsv")
        res <- readRDS(datafile)
        head(res)
        dim(res)

        formatted_res <- res %>%
            mutate(cell_type = ct,
                   summary_type = st, 
                   region_type = rt
            ) %>%
            rename(twas_zscore = zscore,
                   coloc_GLCP = GLCP
            ) 
        head(formatted_res)

        formatted_res
    }
    all_res <- apply(runs, 1, load_intact) %>%
        do.call(rbind, .)
    head(all_res)

    # add gene symbol on to results
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/without_global_pc_cleaning/11_intact/wightman/ensembl_gene_symbols.rds")
    gene_symbols <- readRDS(datafile) %>%
        select(gene=ensembl_gene_id, symbol=hgnc_symbol)
    head(gene_symbols)

    annotated_res <- all_res %>%
        separate(gene, c('gene', 'gene_vers', 'pc'), sep='\\.', fill='right') %>%
        left_join(gene_symbols)
    head(annotated_res)

    # for PCs, select the most significant value for each test
    summarized_res <- annotated_res %>%
        group_by(gene, cell_type, summary_type, region_type, symbol) %>%
        summarize(twas_zscore = max(abs(twas_zscore)),
                  coloc_GLCP = max(coloc_GLCP),
                  intact_pip = max(intact_pip)
        ) 
    head(summarized_res)

    summarized_res
}


plot_data <- function(combined_dat, ad_genes, niagds_genes, region_type){
    # plot the data all together
    head(combined_dat)

    # add columns for significance
    #twas_thresh <- 1.96
    twas_thresh <- 5.45 # genome wide significance (5e-8)
    coloc_thresh <- 0.6
    intact_thresh <- 0.6

    head(combined_dat)
    formatted_dat <- combined_dat %>%
        mutate(twas_sig = abs(twas_zscore) > twas_thresh,
               coloc_sig = coloc_GLCP > coloc_thresh,
               intact_sig = intact_pip > intact_thresh
        ) %>%
        mutate(niagds_gene = symbol %in% niagds_genes$symbol,
               ad_gene = symbol %in% ad_genes$symbol
        ) %>%
        left_join(cell_type_map) %>%
        left_join(summary_type_map)
    head(formatted_dat)

    sig_dat <- formatted_dat %>%
        filter(intact_sig)
    head(sig_dat)
    unique(sig_dat$gene) %>% length() # number of sig genes

    formatted_dat %>%
        arrange(-niagds_gene, -coloc_GLCP) %>%
        head()

    head(combined_dat)
    ranked_dat <- combined_dat %>%
        arrange(-intact_pip) %>%
        #gather('method', 'score', twas_zscore:intact_pip) %>%
        #group_by(method) %>%
        #arrange(-score) %>%
        mutate(rank = row_number()) 
        #select(-score) 
    head(ranked_dat)

    # look for genes that ranked high across all summary type
    #ordered_rank <- ranked_dat %>%
        #group_by(gene, cell_type, summary_type, region_type, symbol) %>%
        #mutate(summed_rank = sum(rank)) %>%
        #spread(method, rank) %>%
        #arrange(summed_rank) 
    #head(ordered_rank)
    #dim(ordered_rank)

    #savefile <- paste0(savedir, "ranked_genes.rds")
    #write_csv(ordered_rank, savefile)

    #filtered_genes <- unique(ranked_dat$gene)
    filtered_genes <- unique(sig_dat$gene)
    filtered_genes
    length(filtered_genes)

    ordered_genes <- data.frame(gene=filtered_genes) %>%
        mutate(gene_rank = row_number())
    
    # get gene position info
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/without_global_pc_cleaning/11_intact/wightman/ensembl_gene_symbols.rds")
    gene_annots <- readRDS(datafile) %>%
        rename(gene = ensembl_gene_id)
    head(gene_annots)

    prioritized_genes <- formatted_dat %>%
        filter(gene %in% filtered_genes) %>%
        left_join(gene_annots) %>%
        left_join(ad_genes)
    head(prioritized_genes)
    dim(prioritized_genes)
    length(unique(prioritized_genes$symbol))

    (savefile <- paste0(savedir, "all_sig_gwas_int_genes.csv"))
    write_csv(prioritized_genes, savefile)

    # check on colocalized genes
    coloc_genes <- prioritized_genes %>%
        filter(coloc_sig)
    head(coloc_genes)

    length(unique(coloc_genes$gene))
    
    stats <- coloc_genes %>%
        filter(region_type == 'full_gene') %>%
        ungroup() %>%
        select(symbol, overallAssociationScore) %>%
        unique() %>%
        filter(!is.na(overallAssociationScore)) %>%
        summarise(min=min(overallAssociationScore),
                  mean=mean(overallAssociationScore),
                  median=median(overallAssociationScore),
                  max=max(overallAssociationScore)
        )
    stats

    head(ad_genes)
    summary(ad_genes$overallAssociationScore)
    qnts <- quantile(ad_genes$overallAssociationScore, probs=seq(0,1,by=0.01))
    (interval <- findInterval(stats$min, qnts, rightmost.closed =TRUE))
    qnts[[interval]]
    names(qnts)[[interval]]


    head(prioritized_genes)
    counts <- prioritized_genes %>%
        ungroup() %>%
        select(symbol, cell_type, summary_type, intact_sig, twas_sig, coloc_sig) %>%
        filter(intact_sig) %>%
        select(symbol, summary_type, intact_sig, cell_type) %>%
        unique() %>%
        count(summary_type, symbol)
    head(counts)

    # looking at counts of sig genes
    table(counts$summary_type, counts$n)
    counts %>%
        filter(summary_type == 'avgs') %>%
        summary()


    head(prioritized_genes)



    #region_type <- "full_gene"
    savefile <- paste0(savedir, region_type, "_gwas_int_results.csv")
    savefile
    write_csv(prioritized_genes, savefile)


    # bring in gwas informatoin
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/02_gwas_prep/02_check_liftover/wightman/formatted_signed_dgap_gwas_scores.txt")
    gwas <- read_tsv(datafile, col_names=c('variant_id', 'gwas_zscore'))
    head(gwas)


    format_table <- prioritized_genes %>%
        select(-twas_zscore, -coloc_GLCP, -intact_pip) %>%
        gather('method', 'sig', twas_sig:intact_sig) %>%
        mutate(method = str_remove(method, "_sig")) %>%
        unite('method_celltype', method, cell_type, sep='_') %>%
        spread(method_celltype, sig) %>%
        left_join(ordered_genes) %>%
        arrange(gene_rank)
    head(format_table)
    dim(format_table)


    # plot manhattan plots
    head(formatted_dat)
    annotated_dat <- formatted_dat %>%
        left_join(gene_annots) 
    head(annotated_dat)

    data_cum <- annotated_dat %>%
        ungroup() %>%
        select(gene, chromosome_name, start_position) %>%
        unique() %>%
        mutate(chromosome_name = as.numeric(chromosome_name)) %>%
        filter(!is.na(chromosome_name)) %>%
        arrange(chromosome_name) %>%
        group_by(chromosome_name) %>%
        mutate(bp=start_position / 1000) %>%
        summarise(max_bp = max(bp, na.rm=TRUE)) %>%
        mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
        select(chromosome_name, bp_add)
    head(data_cum,20)
    tail(data_cum)
    unique(data_cum$chromosome_name)

    ordered_dat <- annotated_dat %>%
        ungroup() %>%
        select(gene, chromosome_name, start_position) %>%
        unique() %>%
        mutate(chromosome_name = as.numeric(chromosome_name)) %>%
        inner_join(data_cum, by='chromosome_name') %>%
        mutate(bp_cum = bp_add + (start_position/1000))
    head(ordered_dat)

    axis_set <- ordered_dat %>%
        group_by(chromosome_name) %>%
        summarize(center = mean(bp_cum))
    head(axis_set)

    head(annotated_dat)
    head(ordered_dat)
    ordered_annots <- annotated_dat %>%
        mutate(chromosome_name = as.numeric(chromosome_name)) %>%
        left_join(ordered_dat)
    head(ordered_annots)

    labeled_annots <- ordered_annots %>%
        mutate(label = ifelse(coloc_GLCP > coloc_thresh, hgnc_symbol, ""))
    head(labeled_annots)

    savefile <- paste0(savedir, region_type, "_colocalization_manhattan.png")
    p <- labeled_annots %>%
        ggplot(aes(x=bp_cum, y=coloc_GLCP, color=as_factor(chromosome_name), size=coloc_GLCP)) +
        geom_hline(yintercept=coloc_thresh, color='grey40', linetype='dashed') +
        geom_point(aes(alpha=0.75)) +
        geom_text_repel(aes(label=label), nudge_x=0.15, nudge_y=0.1, color='black', size=3) +
        facet_grid(cell_type ~ summary_type) +
        scale_x_continuous(label = axis_set$chromosome_name, breaks = axis_set$center) +
        #scale_y_continuous(expand = c(0,0), limist=c(0,ylim)) +
        scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chromosome_name)))) +
        scale_size_continuous(range = c(0.5,3)) +
        theme_bw() +
        theme(strip.background=element_rect(fill="white"),
              legend.position='none',
              panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              #axis.title.y=element_markdown(),
              axis.text.x=element_text(angle=90, size=8, vjust=0.5)
        ) +
        labs(x=NULL, y="Gene-level colocalization probability")
    ggsave(p, file=savefile, width=8, height=7)

    # intact probability
    savefile <- paste0(savedir, region_type, "_intact_manhattan.png")
    head(labeled_annots)
    p <- labeled_annots %>%
        mutate(label = ifelse(intact_pip > intact_thresh, hgnc_symbol, NA)) %>%
        ggplot(aes(x=bp_cum, y=intact_pip, color=as_factor(chromosome_name), size=intact_pip)) +
        geom_hline(yintercept=intact_thresh, color='grey40', linetype='dashed') +
        geom_point(aes(alpha=0.75)) +
        geom_text_repel(aes(label=label), nudge_x=0.15, nudge_y=0.1, color='black', size=3, force=2,
                        force_pull=0.5) +
        facet_grid(cell_type ~ summary_type) +
        scale_x_continuous(label = axis_set$chromosome_name, breaks = axis_set$center, expand=c(0,0)) +
        scale_y_continuous(expand = c(0.01,0.01), limits=c(0,1.1)) +
        #scale_color_manual(values = rep(c("#f56476", "#70d6ff"), unique(length(axis_set$chromosome_name)))) +
        scale_color_manual(values = rep(c("#0B3954", "#F33F54"), unique(length(axis_set$chromosome_name)))) +
        scale_size_continuous(range = c(0.5,3)) +
        theme_bw() +
        theme(strip.background=element_rect(fill="white"),
              legend.position='none',
              panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              #axis.title.y=element_markdown(),
              axis.text.x=element_text(angle=90, size=8, vjust=0.5)
        ) +
        labs(x=NULL, y="Gene-level causality probability")
    ggsave(p, file=savefile, width=9, height=7)


    # plot the three together
    head(ordered_annots)

    # get labels for twas
    twas_labels <- ordered_annots %>%
        ungroup() %>%
        filter(
               #cell_type == 'oligo_opc',
               summary_type == 'pcs') %>%
        select(hgnc_symbol, twas_zscore, chromosome_name) %>%
        unique() %>%
        arrange(-twas_zscore) %>%
        filter(twas_zscore > 5) %>%
        group_by(chromosome_name) %>%
        top_n(n=5, twas_zscore)
    head(twas_labels)
    dim(twas_labels)

    combined_annots <- ordered_annots %>%
        ungroup() %>%
        filter(
               #cell_type == 'oligo_opc',
               #cell_type != 'bulk',
               summary_type == 'pcs') %>%
        select(cell_type:intact_pip, niagds_gene, hgnc_symbol, chromosome_name, start_position, bp_cum) %>%
        gather('method', 'score', twas_zscore:intact_pip) %>%
        separate(method, c('method', 'score_type'), sep='_') %>%
        mutate(thresh = ifelse(method == 'twas', twas_thresh, intact_thresh)) %>%
        mutate(label = ifelse(hgnc_symbol %in% unique(prioritized_genes$hgnc_symbol), hgnc_symbol, "")) %>%
        mutate(label = ifelse(hgnc_symbol %in% unique(twas_labels$hgnc_symbol) & method == 'twas', hgnc_symbol, label)) %>%
        mutate(label = ifelse(score > thresh & method != 'twas', hgnc_symbol, label)) 
    head(combined_annots)


    (savefile <- paste0(savedir, region_type, "_combined_methods_manhattan.png"))
    head(combined_annots)
    
    unique(combined_annots$score_type)

    keep_labels <- combined_annots %>%
        filter(cell_type == 'neuron',
               summary_type == 'pcs',
               region_type == 'full_gene'
               ) %>%
        filter(score > thresh,
               method == 'intact'
        ) %>%
        select(symbol) %>%
        unique()
    head(keep_labels)

    p <- combined_annots %>%
        mutate(method = case_when(method == "coloc" ~ "Colocalization",
                                  method == "intact" ~ "INTACT",
                                  method == "twas" ~ "TWAS")) %>%
        mutate(method = factor(method, levels=c('Colocalization', 'TWAS', 'INTACT'))) %>%
        filter(cell_type == 'neuron',
               summary_type == 'pcs',
               region_type == 'full_gene'
               ) %>%
        mutate(label = ifelse(#(score > thresh & symbol %in% ad_genes$symbol) | 
                              hgnc_symbol %in% keep_labels$symbol, hgnc_symbol, "")) %>%
        mutate(label_color = ifelse(hgnc_symbol %in% keep_labels$symbol, 'blue', 'black')) %>%
        left_join(cell_type_map) %>%
        ggplot(aes(x=bp_cum, y=score, color=as_factor(chromosome_name), size=score)) +
        facet_wrap(method ~ ., ncol=1, scales='free_y') +
        #facet_grid(method ~ formatted_ct, scales='free_y') +
        geom_hline(aes(yintercept=thresh), color='grey40', linetype='dashed') +
        geom_point(aes(alpha=0.75), size=1) +
        geom_text_repel(aes(label=label), color='black', nudge_x=0.15, nudge_y=0.1, size=5, force=2,
                        force_pull=1, max.overlaps=20) +
        scale_x_continuous(label = axis_set$chromosome_name, breaks = axis_set$center, expand=c(0,0)) +
        #scale_y_continuous(expand = c(0.01,0.01), limits=c(0,1.1)) +
        #scale_color_manual(values = rep(c("#f56476", "#70d6ff"), unique(length(axis_set$chromosome_name)))) +
        scale_color_manual(values = rep(c("#0B3954", "#F33F54"), unique(length(axis_set$chromosome_name)))) +
        #scale_size_continuous(range = c(0.5,3)) +
        theme_bw() +
        theme(strip.background=element_rect(fill="white", color='white'),
              legend.position='none',
              panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              text=element_text(size=18),
              #axis.title.y=element_markdown(),
              axis.text.x=element_text(angle=90, size=11, vjust=0.5),
              panel.border=element_rect(size=2)
        ) +
        labs(x="Chromosome", y="Gene-level association score") 
        #ggtitle("oligo_opc")
    savefile
    ggsave(p, file=savefile, width=9, height=10)


    # making a plain manhattan for slides
    (savefile <- paste0(savedir, "plain_manhattan.png"))

    p <- combined_annots %>%
        mutate(method = case_when(method == "coloc" ~ "Colocalization",
                                  method == "intact" ~ "INTACT",
                                  method == "twas" ~ "TWAS")) %>%
        mutate(method = factor(method, levels=c('Colocalization', 'TWAS', 'INTACT'))) %>%
        filter(cell_type == 'neuron',
               summary_type == 'pcs',
               region_type == 'full_gene',
               method == 'TWAS'
               ) %>%
        mutate(label = ifelse(#(score > thresh & symbol %in% ad_genes$symbol) | 
                              hgnc_symbol %in% keep_labels$symbol, hgnc_symbol, "")) %>%
        mutate(label_color = ifelse(hgnc_symbol %in% keep_labels$symbol, 'blue', 'black')) %>%
        left_join(cell_type_map) %>%
        ggplot(aes(x=bp_cum, y=score, color=as_factor(chromosome_name), size=score)) +
        #facet_wrap(method ~ ., ncol=1, scales='free_y') +
        #facet_grid(method ~ formatted_ct, scales='free_y') +
        geom_hline(aes(yintercept=thresh), color='grey40', linetype='dashed') +
        geom_point(aes(alpha=0.75), size=1) +
        #geom_text_repel(aes(label=label), color='black', nudge_x=0.15, nudge_y=0.1, size=5, force=2,
                        #force_pull=1, max.overlaps=20) +
        scale_x_continuous(label = axis_set$chromosome_name, breaks = axis_set$center, expand=c(0,0)) +
        #scale_color_manual(values = rep(c("#f56476", "#70d6ff"), unique(length(axis_set$chromosome_name)))) +
        #scale_color_manual(values = rep(c("#0B3954", "#F33F54"), unique(length(axis_set$chromosome_name)))) +
        scale_color_manual(values = rep(c("#3FA9F5", "#000000"), unique(length(axis_set$chromosome_name)))) +
        #scale_size_continuous(range = c(0.5,3)) +
        theme_bw() +
        theme(strip.background=element_rect(fill="white", color='white'),
              legend.position='none',
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              #text=element_text(size=25),
              axis.text=element_blank(),
              axis.title=element_blank(),
              axis.ticks=element_blank(),
              #axis.title.y=element_markdown(),
              #axis.text.x=element_text(angle=90, size=11, vjust=0.5),
              panel.border=element_rect(size=2)
        ) +
        labs(x="Chromosome", y="Gene-level association score") 
        #ggtitle("oligo_opc")
    savefile
    ggsave(p, file=savefile, width=9, height=4)
}

load_ad_genes <- function(){
    datafile <- paste0("/home/eulalio/deconvolution/data/gene_sets/AD_opentargets.csv")
    ad_genes <- read_csv(datafile)
    head(ad_genes)
    ad_genes
}

main <- function(){
    cell_types <- c('bulk', 'astro', 'endo', 'neuron', 'oligo_opc')
    region_types <- c('full_gene', 'gene_body', 'promoters', 'preTSS')
    summary_types <- c('pcs', 'avgs')

    runs <- expand.grid(cell_type=cell_types,
                        region_type=region_types,
                        summary_type=summary_types
    ) %>%
        mutate(rn = row_number()) %>%
        mutate(runname = paste(summary_type, cell_type, region_type, sep='_'))
    head(runs)


    # load ad opentargets genes
    ad_genes <- load_ad_genes()
    datafile <- paste0("/home/eulalio/deconvolution/data/gene_sets/NIAGDS_AD_genes_list.csv")
    niagds_genes <- read_csv(datafile, col_names='symbol')
    head(niagds_genes)


    # load the data
    combined_dat <- load_data(runs)
    head(combined_dat)
    (savefile <- paste0(savedir, "gwas_integration_scores.csv"))
    write_csv(combined_dat, savefile)

    # plot the data
    #region_type <- "all_regions"
    (region_type <- region_types[1])
    plot_data(combined_dat, ad_genes, niagds_genes, region_type)
    
}
