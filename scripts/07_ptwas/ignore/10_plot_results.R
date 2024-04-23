library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(tidyverse)


source("/home/eulalio/deconvolution/new_rosmap/scripts/shared_functions/02_plot_colors.R")

#proportions_type="without_proportions"
proportions_type="with_proportions"

savedir <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/10_plot_results/",
                  proportions_type, "/"
)
dir.create(savedir)

plotdir <- paste0(savedir, "plots/")
dir.create(plotdir)

load_data <- function(summary_type, region_type, cell_type, runname){
    # load the ptwas estimates data
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/09_ptwas_est/",
                       runname, "/", gwas_type, "_final_ptwas_est.txt")
    res <- read_tsv(datafile)

    head(res)

    res
}

load_ad_genes <- function(){
    datafile <- paste0("/home/eulalio/deconvolution/data/gene_sets/AD_opentargets.csv")
    ad_genes <- read_csv(datafile)
    head(ad_genes)
    ad_genes
}


format_res <- function(res, gene_map, ad_genes){
    # explore the results
    head(res)
    #str(res)

    summary(res$cum_pip)
    summary(res$estimated_eff)
    summary(res$num_instruments)
    summary(res$`I^2`)
    dim(res)

    formatted_res <- res %>%
        as.data.frame() %>%
        rename(gene="#Gene") %>%
        separate(gene, c('gene_nover', 'gene_vers', 'pc'), sep='\\.', remove=FALSE, fill='right')
    head(formatted_res)
    #str(formatted_res)

    # map gene IDs to gene symbols
    savefile <- paste0(savedir, "biomart_genes.rds")
    if (!file.exists(savefile)){
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        sig_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 #filters = "ensembl_gene_id", 
                 #values = unique(formatted_res$gene_nover), 
                 mart = ensembl) 
        saveRDS(sig_symbols, savefile)
    } else{
        print("Loading existing genes symbol file")
        sig_symbols <- readRDS(savefile)
    }
    head(sig_symbols)

    res_symbols <- formatted_res %>%
        rename(ensembl_gene_id=gene_nover) %>%
        left_join(sig_symbols)
    head(res_symbols)
    
    # mark whether these are AD genes or not
    head(ad_genes)

    ad_res <- res_symbols %>%
        rename(symbol=hgnc_symbol) %>%
        left_join(ad_genes[c('symbol', 'overallAssociationScore')]) %>%
        mutate(ad_gene = !is.na(overallAssociationScore))
    head(ad_res)

    # add zscore
    zscore_res <- ad_res %>%
        mutate(zscore = estimated_eff / std_err)
    head(zscore_res)

    zscore_res
}

plot_results <- function(res_df){
    # plot the results
    head(res_df)

    summary(res_df$zscore)

    # get the counts of sig results across cell cell_types
    p <- res_df %>%
        #filter(abs(zscore) > 1.65) %>%
        filter(cum_pip > 0.5) %>%
        mutate(formatted_ct = fct_recode(cell_type,
                                         "Bulk"="bulk",
                                         "Astrocyte"="astro",
                                         "Endothelial cell"="endo",
                                         "Neuron"="neuron",
                                         "Oligo/OPC"="oligo_opc"
            )
        ) %>%
        mutate(formatted_ct = fct_relevel(formatted_ct,
                                          "Bulk", "Astrocyte", "Endothelial cell", "Neuron", "Oligo/OPC")) %>%
        mutate(formatted_st = fct_recode(summary_type,
                                         "Avgs"="avgs",
                                         "PCs"="pcs"
            )) %>%
        ggplot() +
        geom_bar(aes(x=formatted_ct, fill=formatted_st),
                 position='dodge'
        ) +
        scale_fill_manual(values=summary_type_colors, name="Summary type") +
        theme_bw() +
        theme(text=element_text(size=20),
              axis.text.x=element_text(angle=30, hjust=1, vjust=1),
              strip.background=element_rect(fill='white', color='black'),
              #panel.spacing=unit(0, 'lines'),
              panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
        ) +
        xlab("Cell type") +
        ylab("Significant genes") +
        ggtitle(paste(gwas_type, region_type))
    (savefile <- paste0(plotdir, gwas_type, "_", region_type, "_ptwas_est_bar_plot.png"))
    ggsave(p, file=savefile, width=12, height=10)

    
    head(res_df)
    p <- res_df %>%
        #filter(abs(zscore) > 1.65) %>%
        filter(cum_pip > 0.5) %>%
        #arrange(-abs(zscore)) %>%
        arrange(-cum_pip) %>%
        head(30) %>%
        mutate(symbol = ifelse(is.na(symbol), gene_id, symbol)) %>%
        mutate(symbol = ifelse(is.na(pc), symbol, paste(symbol, pc, sep='-'))) %>%
        #mutate(symbol = ifelse(niagds_gene, paste0("*", symbol), symbol)) %>%
        mutate(formatted_ct = fct_recode(cell_type,
                                         "Bulk"="bulk",
                                         "Astrocyte"="astro",
                                         "Endothelial cell"="endo",
                                         "Neuron"="neuron",
                                         "Oligo/OPC"="oligo_opc"
            )
        ) %>%
        mutate(formatted_ct = fct_relevel(formatted_ct,
                                          "Bulk", "Astrocyte", "Endothelial cell", "Neuron", "Oligo/OPC")) %>%
        mutate(formatted_st = fct_recode(summary_type,
                                         "Avgs"="avgs",
                                         "PCs"="pcs"
            )
        ) %>%
        mutate(known_ad_gene = ifelse(is.na(overallAssociationScore), "Novel gene", "AD gene")) %>%
        ggplot(aes(x=formatted_ct, y=fct_reorder(symbol, overallAssociationScore), fill=formatted_ct)) +
        geom_tile(aes(alpha=cum_pip)) +
        geom_text(aes(label=round(cum_pip,2))) +
        facet_grid(known_ad_gene ~ formatted_st, scale='free_y', space='free') +
        scale_fill_manual(values=cell_type_colors, name="Cell type") +
        theme_bw() +
        theme(text=element_text(size=20),
              axis.text.x=element_text(angle=30, hjust=1, vjust=1),
              strip.background=element_rect(fill='white', color='black'),
              #panel.spacing=unit(0, 'lines'),
              panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
        ) +
        xlab("Cell type") +
        ylab("Colocalized genes")
    (savefile <- paste0(plotdir, gwas_type, "_", region_type, "_ptwas_est_tile_plot.png"))
    ggsave(p, file=savefile, width=12, height=10)
}

main <- function(){
    cell_types <- c('bulk', 'astro', 'neuron', 'oligo_opc')
    summary_types <- c('avgs', 'pcs')
    region_types <- c('full_gene')
    gwas_types <- c('wightman')

    gt <- gwas_types[1]
    region_type <- region_types[1]
    gwas_type <- gt

    summary_type <- 'pcs'
    region_type <- 'full_gene'
    cell_type <- 'bulk'

    # load ad opentargets genes
    ad_genes <- load_ad_genes()


    runs <- expand.grid(summary_type=summary_types,
                        cell_type=cell_types,
                        region_type=region_types
    )
    head(runs)

    # load data across all runs
    row <- runs[1,]
    process_row <- function(row){
        summary_type <- row[['summary_type']] %>% as.character()
        region_type <- row[['region_type']] %>% as.character()
        cell_type <- row[['cell_type']] %>% as.character()
        
        # load the data
        runname <- paste(summary_type, cell_type, region_type, sep='_')

        print(paste("Processing", runname))

        res <- load_data(summary_type, region_type, cell_type, runname)

        # check res
        formatted_res <- format_res(res, gene_map, ad_genes)
        head(formatted_res)
        #str(formatted_res)


        # add column for this run
        final_res <- formatted_res %>%
            mutate(summary_type=summary_type,
                   region_type=region_type,
                   cell_type=cell_type
            )

        savefile <- paste0(savedir, runname, "_formatted_ptwas_scan_results.rds")
        savefile
        saveRDS(final_res, savefile)

        final_res
    }
    res <- apply(runs, 1, process_row)
    res_df <- do.call(rbind, res)
    head(res_df)

    savefile <- paste0(savedir, gwas_type, region_type, "_combined_ptwas_est_results.csv")
    write_csv(res_df, savefile)

    plot_results(res_df)
} 

