library(data.table)
library(tidyverse)

# identify variants to exclude based on missing haplotype test results


# Initial setup
datafile <- paste0("../../output/05_qtl_analysis/04_mishap_test/all_chroms_hg38.missing.hap")
thresh <- 1e-9
savedir <- paste0("../../output/05_qtl_analysis/05_filter_mishap/")

# Create output directory if not exists
dir.create(savedir)


## Process the missing haplotype results
# Read in the mishap results
hap <- read_table(datafile)
head(hap)

# filter for significance < 1e-9 based on xQTL study
sig <- hap %>%
    filter(P < thresh)
head(sig)
dim(sig)

# Prepare the snps to exclue
exclude_df <- data.frame(snp=sig$SNP)

# Save the exclude variants list
savefile <- paste0(savedir, "exclude_variants_hg38.txt")
write_tsv(exclude_df, savefile, col_names=FALSE)


# create a file to merge the filesets together
bim_files <- paste0("../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_chr", chroms)
bim_df <- data.frame(bims = bim_files)

savefile <- paste0(savedir, "bim_files_hg38.txt")
write_tsv(bim_df, savefile, col_names=FALSE)


# store the refiltered names too
refiltered_bim_files <- paste0("../../output/05_qtl_analysis/06_merge_bims/rosmap_wgs_hg38_chr", chroms, "_refiltered")
refiltered_bim_df <- data.frame(bims = refiltered_bim_files)

savefile <- paste0(savedir, "bim_files_hg38_refiltered.txt")
write_tsv(refiltered_bim_df, savefile, col_names=FALSE)


# get list of samples to filter out to match methylation data
datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds")
matched_clinical <- readRDS(datafile)

keep_samples <- data.frame(fam=0, specimenID=matched_clinical$wgs.specimenID)
savefile <- paste0(savedir, "keep_samples.txt")
write_tsv(keep_samples, savefile, col_names=FALSE)

# get list of all wgs samples
datafile <- ("/home/eulalio/deconvolution/rosmap/data/ROSMAP_data/Metadata/ROSMAP_assay_wholeGenomeSeq_metadata.csv")
meta <- read_csv(datafile)

# Identify the samples to be removed
remove_samples <- data.frame(specimenID=setdiff(meta$specimenID, matched_clinical$wgs.specimenID)) %>%
    unique()

savefile <- paste0(savedir, "remove_samples.txt")
write_csv(remove_samples, savefile, col_names=FALSE)
