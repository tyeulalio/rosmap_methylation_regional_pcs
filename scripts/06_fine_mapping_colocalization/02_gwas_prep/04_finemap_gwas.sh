#!/usr/bin/env bash

TMPDIR=$1
REGIONNUM=$2

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate colocalization

# Fine-map GWAS using DAP-G and QTL pipeline scripts
GWASTYPE="wightman"

SNP_MAP_FILE="/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/02_check_liftover/${GWASTYPE}/gwas_snp_ld_block_annots.tsv"
ZSCORE_FILE="/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/02_check_liftover/${GWASTYPE}/formatted_signed_dgap_gwas_scores.txt"

OUTPUT_DIR="/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/04_finemap_gwas/${GWASTYPE}"
mkdir -p $OUTPUT_DIR

# plink formatted genotype data
PLINK_FILE="/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/03_create_gwas_vcfs/plink/rosmap_wgs_hg38_all_chroms_no_dups"

THREADS=4

python3 /path/to/qtl_pipeline/02_finemapping/finemap_gwas/01_finemap_gwas_with_dapg.py \
    -m $SNP_MAP_FILE  \
    -r $REGIONNUM \
    -z $ZSCORE_FILE \
    -o $OUTPUT_DIR \
    -p $PLINK_FILE \
    -th $THREADS \
    -t $TMPDIR
