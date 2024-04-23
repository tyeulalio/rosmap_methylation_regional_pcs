library(tidyverse)

# check the ptwas scan results


summary_types <- c('avgs')
region_types <- c('preTSS')
cell_types <- c('bulk', 'astro', 'endo', 'neuron', 'oligo_opc')
proportion_type <- "with_proportions"

runs <- expand.grid(
                    region_type=region_types,
                    cell_type=cell_types,
                    summary_type=summary_types
    ) %>%
    mutate(runname = paste(summary_type, cell_type, region_type, sep='_'))
head(runs)

gwas_type <- "wightman"

runs
run <- runs[3,]
process_run <- function(run){
    runname <- run[['runname']]
    print(paste("Processing", runname))
    
    chrom <- 1
    process_chrom <- function(chrom){

        # summary file
        datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/07_ptwas_scan/",
                           gwas_type, "/",
                           proportion_type, "/",
                           runname, "/chr", chrom, "_ptwas_scan.summary_out.txt")
        datafile
        if (!file.exists(datafile)){
            print(paste("No output for", runname, "chromosome", chrom))
            return(NA)
        }

        datafile
        summary_res <- read_tsv(datafile, col_names=c('chrom', 'pos', 'gene', 'n_snps', 'n_classes', 'top_class', 'top_subclass',
                                                      'min_unadj_pval', 'naive_pval', 'pval', 'info'), skip=1,
                                show_col_types =FALSE
        )
        head(summary_res)

        # stratified file
        datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/07_ptwas_scan/",
                           gwas_type, "/",
                           proportion_type, "/",
                           runname, "/chr", chrom, "_ptwas_scan.stratified_out.txt")
        strat_res <- read_tsv(datafile, col_names=c('chrom', 'pos', 'gene', 'class', 'subclass', 'n_snps', 'stat', 'pval', 'info'),
                              skip=1,
                              show_col_types =FALSE
        )
        head(strat_res)

        strat_res %>%
            arrange(pval)

        head(summary_res)
        head(strat_res)

        dim(summary_res)
        dim(strat_res)

        return(list(strat_res=strat_res,
                    summary_res=summary_res))
    }

    chrom_res <- lapply(1:22, process_chrom)
    head(chrom_res)

    if (all(is.na(chrom_res))){
        print(paste("No output for all of", runname))
        return(1)
    }

    full_res <- chrom_res[!is.na(chrom_res)]

    # combine the data together
    combined_strat <- do.call(rbind, lapply(full_res, function(x) x$strat_res))
    head(combined_strat)

    combined_summary <- do.call(rbind, lapply(full_res, function(x) x$summary_res))
    head(combined_summary)


    nrow(combined_summary)
    table(combined_summary$pval < 0.05)

    check_missing = FALSE
    if (check_missing){
        # load genes from previous step
        datadir <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/01_ptwas_weights/with_proportions/", 
                          runname, "/ptwas_weights/")
        genes <- list.files(datadir) %>%
            str_remove('_ptwas_weights.txt')
        head(genes)
        
        # check which genes are missing
        head(combined_summary)
        missing_genes <- setdiff(genes, unique(combined_summary$gene))
        length(missing_genes)
        head(missing_genes)
    }

    # save to output file
    savedir <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/07_ptwas_scan/", 
                      gwas_type, "/",
                      proportion_type, "/",
                      runname, "/") 
    (savefile <- paste0(savedir, "all_chroms_ptwas_scan.summary_out.txt"))
    write_tsv(combined_summary, savefile)

    (savefile <- paste0(savedir, "all_chroms_ptwas_scan.stratified_out.txt"))
    write_tsv(combined_strat, savefile)
print(paste("Num genes:", length(unique(combined_summary$gene))))
    1
    length(unique(combined_summary$gene))
}
runs
process_run(runs[1,])

apply(runs, 1, process_run)
