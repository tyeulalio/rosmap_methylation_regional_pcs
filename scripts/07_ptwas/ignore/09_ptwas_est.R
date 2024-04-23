##########
### This script reads in DAPG output and then runs it through PTWAS estimation.
### At the end of the script, the individual output files are compiled into one main file.
##########

# Load libraries and downstream parameters
library(tidyverse)

main_savedir <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/09_ptwas_est/")
dir.create(main_savedir)

# for testing
gwas_type <- 'wightman'
summary_type <- 'avgs'
region_type <- 'full_gene'
cell_type <- 'astro'

# get command line arguments
args = commandArgs(trailingOnly=TRUE)
summary_type <- args[1] %>% as.character()
region_type <- args[2] %>% as.character()
cell_type <- args[3] %>% as.character()

# direction for current run
runname <- paste(summary_type, cell_type, region_type, sep='_')
runname
savedir <- paste0(main_savedir, runname, "/")
savedir
dir.create(savedir, showWarnings=FALSE)

load_scan_results <- function(){
    # load in ptwas scan results
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/07_ptwas_scan/",
                       gwas_type, 
                       "/with_proportions/",
                       runname, "/all_chroms_ptwas_scan.summary_out.txt")
    scan_res <- read_table(datafile)
    head(scan_res)

    scan_res
}

# format GWAS scores
format_gwas <- function(){
    # -- match map to the gwas scores file
    (gwas_formatted_file <- paste0(main_savedir, gwas_type, "_formatted_gwas.rds"))
    if (file.exists(gwas_formatted_file)){
        print(paste("Loading existing formatted gwas"))
        formatted_gwas <- readRDS(gwas_formatted_file)
    } else{
        # check dap with gwas names
        gwas_file <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/04_formatted_gwas/wightman_estimated_effect_sizes.tsv")
        gwas <- read_tsv(gwas_file)

        #gwas_map <- readRDS(paste0("/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/02_gwas_prep/02_check_liftover/", gwas_type, "/matched_chrom_pos_mapped_to_variant_name.tsv"))

        # get gwas effect data
        #gwas_file <- read_csv(paste0("/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/02_gwas_prep/02_check_liftover/wightman/gwas_snp_ld_block_annots.tsv"))
        #head(gwas_file)

        # create var to join scores to dap
        #head(gwas)
        #gwas_scores <- gwas %>%
            #unite(var, chrom, pos, ref, alt, sep='_', remove=FALSE)
        #head(gwas_scores)
        #head(gwas_map)

        #length(intersect(gwas_scores$var, gwas_map$var))
        formatted_gwas <- gwas %>%
            select(var, beta, SE)

        #formatted_gwas <- gwas_scores %>%
            #inner_join(gwas_map, by="var") %>%
            #select(var, updated_beta, gwas_SE)
        head(formatted_gwas)
        dim(formatted_gwas)

        gwas_formatted_file
        saveRDS(formatted_gwas, gwas_formatted_file)
    }
    head(formatted_gwas)

    formatted_gwas
}


process_genes <- function(formatted_gwas, sig_scan){ 
    # Begin for loop to process each gene's DAP output
    dap_dir <- paste0("/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/04_dapg/with_proportions/",
                      runname, "/")
    dap_files <- list.files(dap_dir)
    head(dap_files)
    length(dap_files)

    # create a data frame of files
    head(sig_scan)
    genes_df <- data.frame(dap_files=dap_files) %>%
        mutate(gene = str_remove(dap_files, ".dapg"),
               gene = str_replace_all(gene, '-', '.')) %>%
        filter(gene %in% sig_scan$gene) %>%
        mutate(rn = row_number())
    head(genes_df)
    dim(genes_df)

    row <- genes_df[1,]
    process_gene <- function(row){
        #For each dap file, select the proper columns and save as txt file to temp directory
        dap_file <- row[['dap_files']] %>% as.character()
        gene <- row[['gene']] %>% as.character()
        rn <- row[['rn']] %>% as.character()

        print(paste("Processing gene", rn, "out of", nrow(genes_df)))
        
        dap_file <- paste0(dap_dir, dap_file)

        # put files that we need into a temp file
        tmpfile <- paste0("/tmp/", runname, "_tmp.txt")
        cmd <- paste("grep \"((\"", dap_file, " | awk '{ if ($5 != -1) print $2,$5,$3,$6,$7}'  | sort -nk2 > ", tmpfile)
        cmd
        system(cmd)

        #Read in files and name appropriately
        dap_output <- read_table(tmpfile, col_names=c("var", "clusterID", "PIP", "eqtl_effect", "std_error_eqtl"))

        # join the gwas data
        ptwas_input <- dap_output %>%
            left_join(formatted_gwas, by=c('var'))
        head(ptwas_input)

        # store this in the temp files now
        ptwas_input_file <- tmpfile
        write_tsv(ptwas_input, ptwas_input_file, col_names=FALSE)
        savefile <- paste0(savedir, gene,  "_ouput.txt")
        savefile

        # run ptwas est
        pip_thresh = 0.5
        ptwas_cmd = paste("./run_ptwas.sh ",
                           ptwas_input_file,
                            gene,
                           savefile,
                           pip_thresh
        )
        ptwas_cmd
        system(ptwas_cmd)
        }

    res <- apply(genes_df, 1, process_gene)
    1
}

clean_files <- function(){
    # remove the empty files
    clean_cmd <- paste0("find ", savedir, " -type f -empty -delete")
    system(clean_cmd)

    # concat all the txt files and only keep the header/column names from one of the files 
    savefile <- paste0(savedir, gwas_type, "_final_ptwas_est.txt")

    create_concat_file_cmd <- paste("touch", savefile)
    create_concat_file_cmd
    system(create_concat_file_cmd)

    concat_cmd <- paste0("cat ", savedir, "*.txt | awk 'NR == 1 || !/^#/' > ", savefile)
    concat_cmd
    system(concat_cmd)
}

main <- function(){
    # load the ptwas scan results
    scan_res <- load_scan_results()

    # keep only the significant scan results
    head(scan_res)
    sig_scan <- scan_res %>%
        mutate(adj_pval = p.adjust(pval, method='BH'))  %>%
        filter(pval< 0.05)
    table(scan_res$pval < 0.05)
    table(sig_scan$adj_pval < 0.05)

    # get the formatted gwas
    formatted_gwas <- format_gwas()

    process_genes(formatted_gwas, sig_scan)

    clean_files()
}

main()
