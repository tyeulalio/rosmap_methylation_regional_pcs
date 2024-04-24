# check AD genes across QTL analysis
# check proximity of sig genes to APOE
library(GenomicRanges)
library(ggplot2)
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

check_qtls <- function(old_runname, total_genes, ad_genes){ 
    # get sig qtls after qtl mapping
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/05_qtl_analysis/12_map_qtls/with_proportions/",
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


    pval <- get_pval(total_genes, ad_genes, qtl_genes)

    pval
}

check_finemapped <- function(runname, total_genes, ad_genes){
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

    (savefile <- paste0(savedir, runname, "_finemapped_genes.csv"))
    write_csv(finemapped_genes, savefile)

    pval <- get_pval(total_genes, ad_genes, finemapped_genes)
    pval
}

check_colocalization <- function(runname, coloc_thresh, gene_symbols){ 
    # get counts of colocalizations
    datafile <- paste0("/path/to/new_rosmap/output/06_fine_mapping_colocalization/04_fastenloc/01_fastenloc_output/with_proportions/", runname, "/wightman/fastqtl_ld.enloc.gene.out")
    gene_colocalizations <- read_table(datafile)
    head(gene_colocalizations)

    head(gene_symbols)

    # filter for high colocalization probabilities
    sig_colocalizations <- gene_colocalizations %>%
        filter(GLCP > coloc_thresh) %>%
        separate(Gene, c('gene_id', 'pc'), sep='-PC', fill='right') %>%
        separate(gene_id, c('gene_id', 'gene_vers'), sep='-') %>%
        left_join(gene_symbols)
    sig_colocalizations

    sig_genes <- unique(sig_colocalizations$Symbol)
    test <- 'coloc'
    eqtls <- check_eqtls(sig_genes, test) 
        
    return(eqtls)

}

check_ptwas <- function(runname, ptwas_thresh, gene_symbols){ 
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
        separate(gene_id, c('gene_id', 'gene_vers'), sep='\\.') %>%
        left_join(gene_symbols)
    head(sig_ptwas)
    dim(sig_ptwas)
    
    sig_genes <- sig_ptwas$Symbol
    test <- 'ptwas'
    check_eqtls(sig_genes, test)
}

check_intact <- function(runname, intact_thresh, gene_symbols){ 
    # get intact counts
    datafile <- paste0("/path/to/new_rosmap/output/07_ptwas/11_intact/wightman/", 
                       runname, "_all_intact_results.tsv") 
    intact <- readRDS(datafile)
    head(intact)

    sig_intact <- intact %>%
        filter(intact_pip > intact_thresh) %>%
        separate(gene, c('gene_id', 'pc'), sep="\\.PC", fill='right') %>%
        select(gene_id) %>%
        separate(gene_id, c('gene_id', 'gene_vers'), sep='\\.') %>%
        left_join(gene_symbols)
    head(sig_intact)

    sig_genes <- sig_intact$Symbol
    test <- 'intact'
    check_eqtls(sig_genes, test)

}

plot_counts <- function(counts_df){
    head(counts_df)

    # format the counts table nicer
    formatted_df <- counts_df %>%
        gather('count_type', 'count', -summary_type, -cell_type, -region_type) %>%
        #separate(runname, c('summary_type', 'cell_type_region_type'), sep='_', extra='merge') %>%
        mutate(cell_type_region_type = str_replace(cell_type_region_type, 'oligo_opc', 'OligoOPC')) %>%
        separate(cell_type_region_type, c('cell_type', 'region_type'), sep='_', extra='merge')
    head(formatted_df)
    unique(formatted_df$summary_type)
    unique(formatted_df$cell_type)
    unique(formatted_df$region_type)

    # save a wide version of the plot
    formatted_wide <- formatted_df %>%
        spread(count_type, count)
    head(formatted_wide)

    savefile <- paste0(savedir, "formatted_counts_wide.csv")
    write_csv(formatted_wide, savefile)

    plotdir <- paste0(savedir, "plots/")
    dir.create(plotdir, showWarnings = FALSE)

    # plot the counts
    head(formatted_df)
    unique(formatted_df$count_type)
    keep_counts <- c('feature_var_tested', 'feature_var_qtl_sig', 'feature_var_finemapped')
    p <- formatted_df %>%
        mutate(count = as.numeric(count)) %>%
        filter(
               #count < 1e5,
               count_type %in% keep_counts) %>%
        ggplot(aes(x=region_type, y=count, fill=cell_type)) +
        #geom_point(size=4, alpha=0.7) +
        geom_bar(position='dodge', stat='identity') +
        #geom_line() +
        facet_grid(count_type ~ summary_type, scales='free') +
        theme_bw() +
        theme(text = element_text(size=20),
              axis.text.x = element_text(angle=60, hjust=1, vjust=1)) 
    (savefile <- paste0(plotdir, "counts_line_point.png"))
    ggsave(p, file=savefile, width=20, height=10)
}

load_ad_genes <- function(){
    datafile <- paste0("/home/eulalio/deconvolution/data/gene_sets/AD_opentargets.csv")
    ad_genes <- read_csv(datafile)
    head(ad_genes)

    # add gene symbol on to results
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/without_global_pc_cleaning/11_intact/wightman/ensembl_gene_symbols.rds")
    gene_symbols <- readRDS(datafile) %>%
        select(gene=ensembl_gene_id, symbol=hgnc_symbol)
    head(gene_symbols)

    ad_genes <- ad_genes %>%
        left_join(gene_symbols)
    head(ad_genes)
    ad_genes
}

#test_genes <- qtl_genes
get_pval <- function(total_genes, ad_genes, test_genes){
    # compute pval for hypergeometric test
    
    # get the full gene list 
    all_genes <- data.frame(gene=total_genes) %>%
        separate(gene, c('gene', 'gene_vers'), sep='\\.') %>%
        mutate(sig_qtl = gene %in% test_genes$gene) %>%
        mutate(ad_gene = gene %in% ad_genes$gene)
    head(all_genes)

    genes_summarized <- all_genes %>%
        mutate(sig_ad = sig_qtl & ad_gene) %>%
        summarize(total_sig = sum(sig_qtl),
                  total_ad = sum(ad_gene),
                  sig_ad = sum(sig_ad)
        ) %>%
        mutate(total_genes = nrow(all_genes)) %>%
        mutate(total_noad = total_genes - total_ad)
    genes_summarized

    x = genes_summarized[['sig_ad']]
    K = genes_summarized[['total_ad']]
    N = genes_summarized[['total_genes']]
    n = genes_summarized[['total_sig']]

    pval <- phyper(x-1, K, N-K, n,
           lower.tail=FALSE
    )

    data.frame(sig_ad=x,
               total_ad=K,
               total_genes=N,
               total_sig=n,
               pval=pval
    )
}

load_centimorgan_map <- function(){
    # load centimorgan loci map
    datafile <- paste0("/home/eulalio/deconvolution/data/centimorgan_map_hg38_withX.txt.gz")
    centi_map <- read_table(datafile)
    names(centi_map) <- c('chr', 'pos', 'rate', 'cm')
    head(centi_map)

     # regions we want to identify 
    interest_loci <- data.frame(var=c('apoe4', 'apoe2'),
                                chrom=c(19,19),
                                pos=c(44908684,44908822),
                                build=c('hg38', 'hg38')
    )
    head(interest_loci)

    head(centi_map)

    var_row <- interest_loci[1,]
    check_var <- function(var_row){
        print(var_row)
        #var_row[['chrom']]<- as.numeric(var_row[['chrom']])
        #var_row[['pos']]<- as.numeric(var_row[['pos']])
        var_chrom <- var_row[['chrom']] %>% as.numeric()
        var_pos <- var_row[['pos']] %>% as.numeric()

        centi_var <- centi_map %>%
            mutate(rn = row_number()) %>%
            filter(chr == var_chrom) %>%
            mutate(dist = abs(pos - var_pos)) %>%
            arrange(dist)
        head(centi_var)

        # grab the cm entry closest to our variant
        (cm_entry = centi_var[1,])
        
        # find the 1cm before and after this cm entry
        head(centi_map)
        centi_dist <- centi_map %>% 
            mutate(cm_dist = cm - cm_entry[[4]]) %>%
            filter(!abs(cm_dist) < 1) %>%
            mutate(downstream = sign(cm_dist) > 0) %>%
            filter(chr == var_chrom) %>%
            group_by(downstream) %>%
            mutate(dist = abs(pos - var_pos)) %>%
            arrange(abs(cm_dist), dist) %>%
            slice(n=5) %>%
            rename(
                   chrom=chr,
                   cm_pos=pos
            )
        head(centi_dist)


        res <- var_row %>%
            as.data.frame() %>%
            left_join(centi_dist)
        res
    }

    res <- check_var(interest_loci[1,])
    res
        
    (savefile <- paste0(savedir, "interest_loci.txt"))
    write_tsv(res, savefile)

    # format for easier use with genomic ranges
    formatted_res <- res %>%
        select(var, chrom, cm_pos, downstream) %>%
        mutate(end = ifelse(downstream, 'end', 'start')) %>%
        select(-downstream) %>%
        spread(end, cm_pos) %>%
        mutate(chrom = paste0("chr", chrom))
    head(formatted_res)

    formatted_res
}

load_region_annots <- function(region_type){
    # load annotations for region type
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/03_summarised_genes/01_summarised_genes/restimate_proportions/c2_covariates/INT/remove_age/noCleaningSex/annotated_regions/", region_type, "_annotations.rds")
    region_annots <- readRDS(datafile)
    head(region_annots)

    # get one start and end per gene
    sub_annots <- region_annots %>%
        select(gene_id, seqnames, start, end) %>%
        group_by(gene_id, seqnames) %>%
        summarize(start=min(start),
                  end=max(end))
    head(sub_annots)

    sub_annots
}


find_nearby_genes <- function(centi_map, region_annots){
    # find genes near the interest loci
    head(centi_map)
    head(region_annots)

    # create DF for APOE 
    apoe_gr <- makeGRangesFromDataFrame(centi_map, keep.extra.columns =TRUE)
    apoe_gr

    # create dataframe for genes
    gene_gr <- makeGRangesFromDataFrame(region_annots, keep.extra.columns = TRUE)
    gene_gr

    # find overlaps
    overlaps <- findOverlaps(query=apoe_gr, subject=gene_gr) %>%
        as.data.frame()
    head(overlaps)

    # grab the genes that fell within the boundaries
    overlapped_genes <- region_annots[overlaps$subjectHits,]
    head(overlapped_genes)
    dim(overlapped_genes)

    # add gene symbols to these
    datafile <- paste0("/home/eulalio/deconvolution/data/gene_biomart_table.csv")
    symbols <- read_csv(datafile)
    head(symbols)

    overlapped_symbols <- overlapped_genes %>%
        separate(gene_id, c('gene_id', 'gene_vers'), sep='\\.') %>%
        left_join(symbols) %>%
        unite(gene_id, gene_id, gene_vers, sep=".") %>%
        select(-description, -biotype, -ensembl_gene_id_version)
    head(overlapped_symbols)

    overlapped_symbols
}

check_eqtls <- function(sig_genes, test){
    head(sig_genes)    

    # load eqtls
    datafile <- paste0("/path/to/data/qtls/rosmap_eqtls_fixed.csv")
    eqtls <- read_csv(datafile, col_types=cols(.default='c', Beta='d', `Std Error`='d'))
    head(eqtls)

    eqtl_match <- eqtls %>%
        mutate(test_type = test) %>%
        filter(gene %in% sig_genes)
    head(eqtl_match)
    
    if (nrow(eqtl_match) > 0){
        print(paste("**FOUND EQTL GENES***"))
        print(nrow(eqtl_match))
        return(eqtl_match)
    }
    NA
}

main <- function(){
    summary_types <- c('avgs', 'pcs')
    region_types <- c('full_gene', 'promoters', 'preTSS', 'gene_body')
    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')

    runs <- expand.grid(summary_type = summary_types,
                        region_type = region_types,
                        cell_type = cell_types
    ) %>%
        mutate(rn = row_number())
    head(runs)

    coloc_thresh = 0.5
    intact_thresh = 0.5
    ptwas_thresh = 5.45 # genome wide significance

    # get gene symbols
    datafile <- paste0("/path/to/data/full_biomart_table.csv")
    gene_symbols <- read_csv(datafile) %>%
        separate(Ensembl_ID, c('gene_id', 'gene_vers'), sep='\\.') %>%
        select(-gene_vers)
    head(gene_symbols)
        
    row <- runs[1,]
    process_run <- function(row){
        summary_type <- row[['summary_type']]
        region_type <- row[['region_type']]
        cell_type <- row[['cell_type']]
        rn <- row[['rn']]

        (runname <- paste(summary_type, cell_type, region_type, sep='_'))
        print(paste("Processing counts for", runname, rn, "out of", nrow(runs), Sys.time()))

        total_genes <- NA
        ad_genes <- NA

        # load region type annotations
        #region_annots <- load_region_annots(region_type)
        # get genes that are in proximity to apoe
        #proximity_genes <- find_nearby_genes(centi_map, region_annots)
        #head(proximity_genes)

        proximity_genes=NA

        # get counts of colocalizations
        print("Getting colocalization counts")
        (colocalization_genes <- check_colocalization(runname, coloc_thresh, gene_symbols))

        # get counts of PTWAS
        print("Getting PTWAS counts")
        (ptwas_genes <- check_ptwas(runname, ptwas_thresh, gene_symbols))

        # get counts of INTACT
        print("Getting INTACT counts")
        (intact_genes <- check_intact(runname, intact_thresh, gene_symbols))

        combined_genes <- do.call(rbind, list(colocalization_genes, ptwas_genes, intact_genes)) 
        combined_genes

        if (all(is.na(combined_genes))){
            print('no results')
            return(NA)
        }


        combined_genes <- combined_genes %>%
            mutate(summary_type = summary_type,
                   cell_type=cell_type,
                   region_type=region_type
            )
        combined_genes
    }

    res <- apply(runs, 1, process_run)
    res_df <- do.call(rbind, res) %>%
        filter(!is.na(gene))
    head(res_df)

    unique(res_df$`Cell type`)
    unique(res_df$cell_type)
    formatted_res <- res_df %>%
        mutate(eqtl_celltype = fct_recode(`Cell type`,
                                         micro='Microglia',
                                         neuron='Neurons',
                                         endo='Endothelial cells',
                                         astro='Astrocytes',
                                         oligo_opc='Oligodendroglia'
            )
        )
    head(formatted_res)
    
    matched_res <- formatted_res %>%
        filter(eqtl_celltype == cell_type)
    dim(matched_res)
    matched_res

    (savefile <- paste0(savedir, "matched_celltype_eqtl_gwas_integration_results.csv"))
    write_csv(matched_res, savefile)

    1
}

main()

