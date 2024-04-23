#!/usr/bin/env bash

GWASTYPE=$1

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate colocalization

# format the fine=mapped GWAS output properly as FastEnloc input for colocalization
OUTPUT_DIR="/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/06_formatted_finemapped_gwas/${GWASTYPE}"
mkdir -p $OUTPUT_DIR

# fine-mapped dap directory
DAP_DIR="/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/04_finemap_gwas/${GWASTYPE}/gwas_dapg"

python3 /path/to/qtl_pipeline/02_finemapping/finemap_gwas/02_create_fastqtl_gwas_input.py \
    -o $OUTPUT_DIR \
    -d $DAP_DIR
