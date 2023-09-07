#!/bin/bash

# Pipeline for running PLINK and PLINK2 get necessary files for downstream steps
# This script is designed to run on SCG

# Load required modules
module load plink
module load plink2

# define common directories and files
OUTPUT_DIR="../../output/05_qtl_analysis"
DATA_DIR="../../data/rosmap_wgs_harmonization_plink"



# Step 1: Identify multiallelic sites
# try to merge files, this will throw an error because of multiallelic sites
# use the output of this to remove those sites
echo "-- Finding multiallelic sites --"
plink \
    --make-bed \
    --biallelic-only \
    --merge-list "${OUTPUT_DIR}/05_filter_mishap/bim_files_hg38.txt" \
    --exclude "${OUTPUT_DIR}/05_filter_mishap/exclude_variants_hg38.txt" \
    --keep "${OUTPUT_DIR}/05_filter_mishap/keep_samples.txt" \
    --out "${OUTPUT_DIR}/06_merge_bims/all_chroms_hg38"

# Step 2: Remove multiallelic sites from individual chromosome files
for CHROM in {1..22}
do
    echo "removing multiallelic sites for chromosome ${CHROM}"
        plink2 \
            --pfile "${DATA_DIR}/rosmap_wgs_hg38_chr${CHROM}" \
            --make-pgen \
            --exclude "${OUTPUT_DIR}/06_merge_bims/all_chroms_hg38-merge.missnp" \
            --out "${OUTPUT_DIR}/06_merge_bims/rosmap_wgs_hg38_chr${CHROM}_refiltered"
done

# Step 3: Merge the cleaned chromosome files
# Make both p-gen and bed files
echo "-- Merging files --"
plink2 \
    --make-pgen \
    --pfile "${DATA_DIR}/rosmap_wgs_hg38_all_chroms" \
    --exclude "${OUTPUT_DIR}/05_filter_mishap/exclude_variants_hg38.txt" \
    --keep "${OUTPUT_DIR}/05_filter_mishap/keep_samples.txt" \
    --out "${OUTPUT_DIR}/06_merge_bims/all_chroms_hg38"

plink2 \
    --make-bed \
    --pfile "${DATA_DIR}/rosmap_wgs_hg38_all_chroms" \
    --exclude "${OUTPUT_DIR}/05_filter_mishap/exclude_variants_hg38.txt" \
    --keep "${OUTPUT_DIR}/05_filter_mishap/keep_samples.txt" \
    --out "${OUTPUT_DIR}/06_merge_bims/all_chroms_hg38"

# Step 4: Create counts of reference alleles
echo "-- Creating counts file"
plink2 \
    --pfile "${OUTPUT_DIR}/06_merge_bims/${CHROMDIR}all_chroms_hg38" \
    --out "${OUTPUT_DIR}/06_merge_bims/${CHROMDIR}all_chroms_ref_counts_hg38" \
    --export A-transpose
    
# Step 5: Create a PED file for the merged dataset
echo "-- Creating ped file"
plink2 \
    --pfile "${OUTPUT_DIR}/06_merge_bims/${CHROMDIR}all_chroms_hg38" \
    --out "${OUTPUT_DIR}/06_merge_bims/${CHROMDIR}all_chroms_ref_counts_hg38" \
    --export ped
