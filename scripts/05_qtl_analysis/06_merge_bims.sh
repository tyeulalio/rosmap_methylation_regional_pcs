#!/bin/bash

# using plink on scg
module load plink
module load plink2

# run plink to generate files 
# need to run all three parts here, some are commented out during testing

# select chrom18_only is we are just testing
CHROMDIR=""

# now merge the files
echo "-- Merging files --"
plink2 \
    --make-pgen \
    --pfile "../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms" \
    --exclude "../../output/05_qtl_analysis/05_filter_mishap/${CHROMDIR}exclude_variants_hg38.txt" \
    --keep "../../output/05_qtl_analysis/05_filter_mishap/keep_samples.txt" \
    --out "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_hg38"

plink2 \
    --make-bed \
    --pfile "../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms" \
    --exclude "../../output/05_qtl_analysis/05_filter_mishap/${CHROMDIR}exclude_variants_hg38.txt" \
    --keep "../../output/05_qtl_analysis/05_filter_mishap/keep_samples.txt" \
    --out "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_hg38"

# create counts of reference allele
echo "-- Creating counts file"
plink2 \
    --pfile "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_hg38" \
    --out "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_ref_counts_hg38" \
    --export A-transpose

#create ped file
echo "-- Creating ped file"
plink2 \
    --pfile "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_hg38" \
    --out "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_ref_counts_hg38" \
    --export ped
