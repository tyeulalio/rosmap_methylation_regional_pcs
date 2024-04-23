library(INTACT)
library(biomaRt)
library(RColorBrewer)
library(tidyverse)

# run intact on the ptwas and colocalization output

gwas_type <- 'wightman'

savedir <- paste0("/path/to/output/")
    

load_ptwas <- function(runname, gwas_type){
    # load ptwas data
    datafile <- paste0("/path/to/new_rosmap/output/07_ptwas/07_ptwas_scan/", gwas_type, "/with_proportions/",
                       runname, "/all_chroms_ptwas_scan.stratified_out.txt")
    ptwas <- read_table(datafile)
    head(ptwas)
    ptwas
}

load_fastenloc <- function(runname, gwas_type){
    # load fastenloc GLCP data
    datafile <- paste0("/path/to/new_rosmap/output/06_fine_mapping_colocalization/04_fastenloc/01_fastenloc_output/with_proportions/",
                       runname, "/", gwas_type, "/fastqtl_ld.enloc.gene.out")
    fastenloc <- read_table(datafile)
    head(fastenloc)
    fastenloc
}

run_intact <- function(ptwas, fastenloc, runname){
    # run intact
    head(ptwas) 
    head(fastenloc)

    summary(ptwas$stat)

    # join the columns that we want
    sub_ptwas <- ptwas %>%
        select(gene, zscore=stat)
    head(sub_ptwas)

    # match the gene names 
    sub_fastenloc <- fastenloc %>%
        select(gene=Gene, GLCP) %>%
        mutate(gene = str_replace_all(gene, "-", "\\."))
    head(sub_fastenloc)

    ptwas_fastenloc <- sub_ptwas %>%
        inner_join(sub_fastenloc)
    head(ptwas_fastenloc)
    dim(ptwas_fastenloc)

    res <- intact(GLCP_vec=ptwas_fastenloc$GLCP,
                  z_vec=ptwas_fastenloc$zscore,
                  prior_fun=linear
    )
    head(res)
    summary(res)

    # combine results
    intact_res <- cbind(ptwas_fastenloc, intact_pip=res) %>%
        arrange(-intact_pip)
    head(intact_res)

    savefile <- paste0(savedir, runname, "_all_intact_results.tsv")
    saveRDS(intact_res, savefile)

    # filter for high PIPs
    thresh <- 0.5
    sig_res <- intact_res %>%
        filter(intact_pip > thresh) %>%
        separate(gene, c('gene_novers', 'vers', 'pc'), sep='\\.', fill='right', remove=FALSE)
    head(sig_res)

    if (nrow(sig_res) == 0) return(NA)

    # attach gene names
    savefile <- paste0(savedir, "ensembl_gene_symbols.rds")
    if (!file.exists(savefile)){
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        all_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name",
                                            "start_position", "end_position"), 
                 #filters = "ensembl_gene_id", 
                 #values = unique(sig_res$gene_novers), 
                 mart = ensembl) 
        head(all_symbols)
        saveRDS(all_symbols, savefile)
    } else{
        print("Loading exising gene symbols file")
        all_symbols <- readRDS(savefile)
    }
    head(all_symbols)

    # add the symbols to the results
    formatted_res <- sig_res %>%
        rename(ensembl_gene_id=gene_novers) %>%
        left_join(all_symbols) %>%
        select(-ensembl_gene_id, -vers)
    head(formatted_res)
    dim(formatted_res)
    num_sig <- length(unique(formatted_res$hgnc_symbol))

    print(paste("Found", num_sig, "genes"))
    print(head(formatted_res))

    formatted_res
}

load_ad_genes <- function(){
    datafile <- paste0("/path/to/data/gene_sets/AD_opentargets.csv")
    ad_genes <- read_csv(datafile)
    head(ad_genes)
    ad_genes
}


plot_results <- function(res_df, region_type){
    # plot results

    # load ad opentargets genes
    ad_genes <- load_ad_genes()

    # connect ad_Gene annotations to the results
    datafile <- paste0("/path/to/data/gene_sets/NIAGDS_AD_genes_list.csv")
    niagds_genes <- read_csv(datafile, col_names='symbol')
    head(niagds_genes)


    # format the data for plotting
    head(res_df)
    formatted_gene_res <- res_df %>%
        rename(symbol=hgnc_symbol) %>%
        left_join(ad_genes) %>%
        select(gene, pc:symbol, overallAssociationScore, summary_type:region_type) %>%
        mutate(niagds_gene = symbol %in% niagds_genes$symbol)
    head(formatted_gene_res)

    thresh=0.6
    p <- formatted_gene_res %>%
        filter(intact_pip > thresh) %>%
        mutate(symbol = ifelse(is.na(symbol), gene, symbol)) %>%
        mutate(symbol = ifelse(is.na(pc), symbol, paste(symbol, pc, sep='-'))) %>%
        mutate(symbol = ifelse(niagds_gene, paste0("*", symbol), symbol)) %>%
        mutate(known_ad_gene = ifelse(is.na(overallAssociationScore), "Novel gene", "AD gene")) %>%
        left_join(cell_type_map) %>%
        left_join(region_type_map) %>%
        left_join(summary_type_map) %>%
        ggplot(aes(x=formatted_ct, y=fct_reorder(symbol, niagds_gene), fill=formatted_ct)) +
        geom_tile(aes(alpha=intact_pip)) +
        geom_text(aes(label=round(intact_pip,2))) +
        facet_grid(known_ad_gene ~ formatted_st, scale='free_y', space='free_y') +
        scale_fill_manual(values=cell_type_colors, name="Cell type") +
        guides(alpha=guide_legend("INTACT PIP")) +
        theme_bw() +
        theme(text=element_text(size=20),
              axis.text.x=element_text(angle=30, hjust=1, vjust=1),
              strip.background=element_rect(fill='white', color='black'),
              #panel.spacing=unit(0, 'lines'),
              panel.border=element_rect(color='black', fill=NA, linewidth=0.5),
        ) +
        xlab("Cell type") +
        ylab("Putative causal genes") +
        ggtitle(paste(region_type, "intact PIP >", thresh))
    (savefile <- paste0(plot_dir, gwas_type, "_", region_type, "_thresh", thresh, "_intact_pip_tile_plot.png"))
    ggsave(p, file=savefile, width=12, height=10)
}

main <- function(){
    summary_types <- c('avgs', 'pcs')
    region_types <- c('preTSS')
    cell_types <- c('bulk', 'astro', 'endo', 'neuron', 'oligo_opc')


    runs <- expand.grid(summary_type=summary_types,
                        region_type=region_types,
                        cell_type=cell_types
    )
    runs

    row <- runs[1,]
    process_run <- function(row){
        summary_type <- row[['summary_type']] %>% as.character()
        cell_type <- row[['cell_type']] %>% as.character()
        region_type <- row[['region_type']] %>% as.character()

        runname <- paste(summary_type, cell_type, region_type, sep='_')

        savefile <- paste0(savedir, runname, "_formatted_sig_intact_results.tsv")

        print(paste("Processing", runname))

        # load ptwas data
        ptwas <- load_ptwas(runname, gwas_type)

        # load fastenloc data
        fastenloc <- load_fastenloc(runname, gwas_type)

        res <- run_intact(ptwas, fastenloc, runname)

        if (all(is.na(res))) return(NA)

        savefile <- paste0(savedir, runname, "_formatted_sig_intact_results.tsv")
        savefile
        write_tsv(res, savefile)
    

        res %>%
            mutate(summary_type=summary_type,
                   cell_type=cell_type,
                   region_type=region_type
            )
    }

    res <- apply(runs, 1, process_run)
    res_df <- do.call(rbind, res) %>%
        filter(!is.na(summary_type))
    head(res_df)

    dim(res_df)
    
    region_type <- region_types[1]
    plot_results(res_df, region_type)

}
