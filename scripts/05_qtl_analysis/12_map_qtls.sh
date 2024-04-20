#!/usr/bin/env bash

# map qtls using qtl pipeline script

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate tensorqtl

SUMMARY_TYPE=$1
REGION_TYPE=$2
CELL_TYPE=$3
RUNNAME="${CELL_TYPE}_${REGION_TYPE}_${SUMMARY_TYPE}"

#PROPORTION_TYPE="without_proportions"
PROPORTION_TYPE="with_proportions"

PHENO_FILE="${HOME}/deconvolution/new_rosmap/output/05_qtl_analysis/11_formatted_qtl_input/${PROPORTION_TYPE}/${RUNNAME}_matched_pheno.bed.gz"
COVS_FILE="${HOME}/deconvolution/new_rosmap/output/05_qtl_analysis/11_formatted_qtl_input/${PROPORTION_TYPE}/${RUNNAME}_matched_covariates.tsv"
PLINK_PATH="${HOME}/deconvolution/new_rosmap/output/05_qtl_analysis/09_remove_outliers/all_chroms_hg38"

OUTPUT_DIR="${HOME}/deconvolution/new_rosmap/output/05_qtl_analysis/12_map_qtls/${PROPORTION_TYPE}/${RUNNAME}"
mkdir -p $OUTPUT_DIR

# select the cis-window
CIS_WINDOW=1000000
# adjust cis window for cpgs
if [ $SUMMARY_TYPE == 'cpgs' ]
then
    CIS_WINDOW=500000
fi


python3 ${HOME}/qtl_pipeline/01_qtl_mapping/01_run_tensorqtl.py \
    -output $OUTPUT_DIR \
    -pheno $PHENO_FILE \
    -covariates $COVS_FILE \
    -plink $PLINK_PATH \
    -window $CIS_WINDOW
