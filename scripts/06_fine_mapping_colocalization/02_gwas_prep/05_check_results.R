library(tidyverse)

# check for missing finemapped files

gwas_type <- "wightman"
datadir <- paste0("/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/04_finemap_gwas/", 
                  gwas_type, "/gwas_dapg/")
dapg_files <- list.files(datadir)
length(dapg_files)
head(dapg_files)

head(dapg_files)

ld_blocks <- str_remove(dapg_files, "ld_") %>%
    str_remove("\\.dapg")
head(ld_blocks)

setdiff(1:1361, ld_blocks)
