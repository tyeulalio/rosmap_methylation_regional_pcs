#!/bin/bash

# using plink on scg
module load plink2

# use chrom18_only for testing
CHROMDIR=""
#CHROMDIR="chrom18_only/"

# run plink to generate files 

plink2 \
    --make-bed \
    --output-chr chrM \
    --pfile "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_hg38" \
    --remove "../../output/05_qtl_analysis/08_format_eigenstrat/${CHROMDIR}population_outliers_hg38.txt" \
    --out "../../output/05_qtl_analysis/09_remove_outliers/${CHROMDIR}all_chroms_hg38" 
    
