#!/bin/bash

# liftover the VCFs from hg37 to hg38
# using GATK liftoverVCF

# make sure to load the GATK module
module load gatk/4.0.10.0

DATAFILE="/home/eulalio/deconvolution/data/ad_gwas/wightman_etal_2021_results.vcf"

gatk LiftoverVcf \
    --INPUT="${DATAFILE}" \
    --OUTPUT="/path/to/output/01_liftover_gwas/wightman_etal_2021_results_hg38.vcf" \
    --VERBOSITY=DEBUG \
    --CHAIN="/path/to/liftOver/chains/hg19ToHg38.over.chain.gz" \
    --REFERENCE_SEQUENCE="/reference/RefGenomes/GATK_Resource_Bundle/hg38/hg38.fa" \
    --REJECT="/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/01_liftover_gwas/rejected_records.vcf" \
    --MAX_RECORDS_IN_RAM 10000 \
    --RECOVER_SWAPPED_REF_ALT
