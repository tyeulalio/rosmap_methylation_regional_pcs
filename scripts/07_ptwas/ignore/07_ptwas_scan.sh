#!/usr/bin/env bash

# run PTWAS scan

GWAS_FILE=$1
PTWAS_WEIGHT_FILE=$2
LD_PANEL_FILES=$3
PREFIX=$4


echo "GWAS input file: ${GWAS_FILE}"
echo "PTWAS weight input file: ${PTWAS_WEIGHT_FILE}"
echo "LD PANEL input file: ${LD_PANEL_FILES}"

echo "Writing output to ${OUTPUT_FILE}"

time GAMBIT \
    --gwas $GWAS_FILE \
    --betas $PTWAS_WEIGHT_FILE \
    --ldref $LD_PANEL_FILES \
    --ldref-only \
    --prefix $PREFIX

