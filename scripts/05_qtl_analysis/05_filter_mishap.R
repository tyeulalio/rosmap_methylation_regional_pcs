library(data.table)
library(tidyverse)

# identify variants to exclude based on missing haplotype test results


datafile <- paste0("../../output/05_qtl_analysis/04_mishap_test/all_chroms_hg38.missing.hap")
hap <- read_table(datafile)
head(hap)

# filter for significance < 1e-9 based on xQTL study
thresh <- 1e-9
sig <- hap %>%
    filter(P < thresh)
head(sig)
dim(sig)

# save the snps to exclue
exclude <- data.frame(snp=sig$SNP)
head(exclude)

exclude_df <- exclude

head(exclude_df)
dim(exclude_df)


savedir <- paste0("../../output/05_qtl_analysis/05_filter_mishap/")
dir.create(savedir, showWarnings=FALSE)

savefile <- paste0(savedir, "exclude_variants_hg38.txt")
savefile
write_tsv(exclude_df, savefile, col_names=FALSE)


# create a file to merge the filesets together
bim_files <- paste0("../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_chr", chroms)
head(bim_files)

bim_df <- data.frame(bims = bim_files)
head(bim_df)

savefile <- paste0(savedir, "bim_files_hg38.txt")
write_tsv(bim_df, savefile, col_names=FALSE)


# store the refiltered names too
bim_files <- paste0("../../output/05_qtl_analysis/06_merge_bims/rosmap_wgs_hg38_chr", chroms, "_refiltered")
head(bim_files)

bim_df <- data.frame(bims = bim_files)
head(bim_df)

savefile <- paste0(savedir, "bim_files_hg38_refiltered.txt")
write_tsv(bim_df, savefile, col_names=FALSE)


# get list of samples to filter out to match methylation data
datafile <- paste0("/path/to/output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds")
matched_clinical <- readRDS(datafile)
head(matched_clinical)

keep_samples <- data.frame(fam=0, specimenID=matched_clinical$wgs.specimenID)
head(keep_samples)
dim(keep_samples)
savefile <- paste0(savedir, "keep_samples.txt")
write_tsv(keep_samples, savefile, col_names=FALSE)

# get list of all wgs samples
datafile <- ("/home/eulalio/deconvolution/rosmap/data/ROSMAP_data/Metadata/ROSMAP_assay_wholeGenomeSeq_metadata.csv")
meta <- read_csv(datafile)
head(meta)

remove_samples <- data.frame(specimenID=setdiff(meta$specimenID, matched_clinical$wgs.specimenID)) %>%
    unique()
head(remove_samples)
nrow(remove_samples)
nrow(meta)
nrow(matched_clinical)
length(unique(remove_samples$specimenID))

savefile <- paste0(savedir, "remove_samples.txt")
savefile
write_csv(remove_samples, savefile, col_names=FALSE)

1
