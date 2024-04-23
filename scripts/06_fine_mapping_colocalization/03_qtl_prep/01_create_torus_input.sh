#!/usr/bin/env bash

# create torus input from TensorQTL output

# get command line arguments
SUMMARY_TYPE=$1
CELL_TYPE=$2
REGION_TYPE=$3
PROPORTION_TYPE=$4

# output directory
OUTPUT_DIR="${HOME}/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/01_torus_input/${PROPORTION_TYPE}"
mkdir -p $OUTPUT_DIR

RUNNAME="${SUMMARY_TYPE}_${CELL_TYPE}_${REGION_TYPE}"

# input qtl files
QTL_PREFIX="${CELL_TYPE}_${REGION_TYPE}_${SUMMARY_TYPE}"
set -f #disable star expansion
QTL_FILES="${HOME}/deconvolution/new_rosmap/output/05_qtl_analysis/12_map_qtls/${PROPORTION_TYPE}/${QTL_PREFIX}/cis_qtls.cis_qtl_pairs.chr*.parquet"


python3 ${HOME}/qtl_pipeline/02_finemapping/finemap_qtls/01_create_torus_input.py \
    -o $OUTPUT_DIR \
    -p $RUNNAME \
    -qtl $QTL_FILES
