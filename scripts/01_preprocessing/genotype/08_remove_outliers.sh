#!/bin/bash

# Script to run PLINK2 for removing population outliers in genetic data analysis

# Load the PLINK2 module (only applicable if you're using a module-based system like Lmod)
module load plink2

# Directory variables for better code readability
INPUT_DIR="../../output/05_qtl_analysis/06_merge_bims/"
OUTPUT_DIR="../../output/05_qtl_analysis/09_remove_outliers/"
OUTLIER_FILE="../../output/05_qtl_analysis/08_format_eigenstrat/population_outliers_hg38.txt"

# Input and output filenames
INPUT_PFILE="all_chroms_hg38"
OUTPUT_FILE="all_chroms_hg38"

# Run PLINK2 to remove outliers and create new bed files
plink2 \
    --make-bed \                              # Create PLINK binary bed file
    --output-chr chrM \                       # Output chromosome names as 'chrM' for mitochondrial
    --pfile "${INPUT_DIR}${INPUT_PFILE}" \    # Specify input .pfile
    --remove "$OUTLIER_FILE" \                # List of samples to remove
    --out "${OUTPUT_DIR}${OUTPUT_FILE}"       # Specify output file name
