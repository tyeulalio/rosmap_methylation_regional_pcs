#!/bin/bash

# using plink on scg
module load plink

# run plink to generate files for 

echo "Processing chromosome ${CHROM}"
plink \
    --bfile "/home/eulalio/deconvolution/new_rosmap/data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms" \
    --out "../../output/05_qtl_analysis/04_mishap_test/all_chroms_hg38" \
    --test-mishap 
