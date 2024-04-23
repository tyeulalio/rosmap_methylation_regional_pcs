library(tidyverse)

# check results 

# check 01 results, some are missing
# get full list of genes that are being run
dapdir <- paste0("/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/04_dapg/avgs_bulk_full_gene/")
dapfiles <- list.files(dapdir)
head(dapfiles)

genes <- str_remove(dapfiles, ".dapg")
head(genes)

# get ptwas output files
ptwas_dir <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/01_ptwas_weights/")
ptwas_files <- list.files(ptwas_dir)
head(ptwas_files)

# remove file name and replace period 
ptwas_genes <- str_remove(ptwas_files, "_ptwas_weights.txt") %>%
    str_replace("\\.", "-")
head(ptwas_genes)

length(intersect(ptwas_genes, genes))
length(genes)
length(ptwas_genes)

missing_genes <- setdiff(genes, ptwas_genes)
length(missing_genes)
head(missing_genes)
