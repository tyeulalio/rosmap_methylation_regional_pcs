#!/bin/bash

# PLINK Genotype Analysis Pipeline for Research Paper

# Description:
# This script processes genomic data using PLINK, specifically focusing on 
# the dataset 'rosmap_wgs_hg38_all_chroms'. It's tailored for research reproducibility.

# Usage:
# ./script_name.sh 

# Prerequisites:
# - Ensure PLINK is installed and accessible through module system or globally.
# - Ensure 'rosmap_wgs_hg38_all_chroms' dataset is present at the specified path.

# Exit on error
set -e

# Check for arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <chromosome_number>"
    exit 1
fi

# Parameters
CHROM="$1"
INPUT_DATASET="/home/eulalio/deconvolution/new_rosmap/data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms"
OUTPUT_DIR="../../output/05_qtl_analysis/04_mishap_test/all_chroms_hg38"

# Load necessary tools
echo "Loading required modules..."
module load plink

# PLINK Processing
echo "Processing chromosome ${CHROM}..."
plink \
    --bfile "${INPUT_DATASET}" \
    --out "${OUTPUT_DIR}" \
    --test-mishap 

echo "Processing completed for chromosome ${CHROM}."

# End of Script

