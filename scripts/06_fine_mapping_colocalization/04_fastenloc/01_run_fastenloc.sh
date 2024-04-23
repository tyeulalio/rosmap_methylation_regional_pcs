#!/usr/bin/env bash

# run fastenloc

SUMMARYTYPE=$1
CELLTYPE=$2
REGIONTYPE=$3
GWASTYPE=$4
NUM_GWAS_VARS=$5


RUNNAME="${SUMMARYTYPE}_${CELLTYPE}_${REGIONTYPE}"


PROPORTION_TYPE="with_proportions"

OUTDIR="${HOME}/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/04_fastenloc/01_fastenloc_output/${PROPORTION_TYPE}/${RUNNAME}/${GWASTYPE}"
mkdir -p $OUTDIR

# -total_variants = number of GWAS variants tested

QTL_FILE="${HOME}/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/05_fastqtl_qtl_annotations/${PROPORTION_TYPE}/${RUNNAME}/fastenloc.qtl.annotation.vcf.gz" 
GWAS_FILE="${HOME}/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/02_gwas_prep/06_formatted_finemapped_gwas/${GWASTYPE}/fastqtl_gwas_input_pips.vcf.gz"

fastenloc \
    -total_variants $NUM_GWAS_VARS \
    -thread 1 \
    -prefix "${OUTDIR}/fastqtl_ld" \
    -eqtl $QTL_FILE \
    -gwas $GWAS_FILE
