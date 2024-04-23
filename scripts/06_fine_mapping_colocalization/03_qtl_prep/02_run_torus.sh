#!/usr/bin/env bash

# run torus to get SNP PIPs


# get command line arguments
SUMMARY_TYPE=$1
CELL_TYPE=$2
REGION_TYPE=$3
PROPORTION_TYPE=$4

# create the main output directory
OUTPUT_DIR="${HOME}/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/02_torus"
mkdir -p $OUTPUT_DIR

RUNNAME="${SUMMARY_TYPE}_${CELL_TYPE}_${REGION_TYPE}"

# create the output sub-directory
OUTPUT_DIR="${OUTPUT_DIR}/${PROPORTION_TYPE}"
mkdir -p $OUTPUT_DIR
OUTPUT_DIR="${OUTPUT_DIR}/${RUNNAME}"
mkdir -p $OUTPUT_DIR


INPUT_FILE="${HOME}/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/01_torus_input/${PROPORTION_TYPE}/${RUNNAME}_fastqtl_single_snp_output.tsv.gz"


# call the run_torus script from the qtl pipeline
${HOME}/qtl_pipeline/02_finemapping/finemap_qtls/02_run_torus.sh \
    $INPUT_FILE \
    $OUTPUT_DIR
