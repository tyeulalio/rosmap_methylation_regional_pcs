#!/bin/bash

# Liftover Genotype VCF Pipeline for Research Paper

# Description:
# This script lifts over genotype VCFs from hg17 to hg19 using the GATK LiftoverVCF tool. 
# It concatenates, normalizes, and performs liftover on the dataset.

# Usage:
# ./script_name.sh

# Prerequisites:
# - Ensure GATK and bcftools are installed and accessible through the module system or globally.
# - Necessary dataset files and reference genomes should be available at the specified paths.

# Exit on error
set -e

# Parameters
ROSMAP_VCF_PATH="/path/to/rosmap_wgs_harmonization_variant_calling"
HEADER_PATH="/home/eulalio/deconvolution/data/reference/vcf_header_chr_lines.txt"
FASTA_REF_HG19="/reference/RefGenomes/GATK_Resource_Bundle/hg19/hg19.fa" 
OUTPUT_FILE="../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_hg38_all_chroms.vcf" 
CHAIN_FILE="/home/eulalio/shared/liftOver/chains/hg19ToHg38.over.chain.gz" 
FASTA_REF_HG38="/reference/RefGenomes/GATK_Resource_Bundle/hg38/hg38.fa" 
REJECTED_FILE="../../output/05_qtl_analysis/01_liftover_vcf/rejected_records_all_chroms.vcf" 
NORMALIZED_OUTPUT="../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_all_chroms_normalized.vcf" 

# Temporary storage for operations
SAVETMP="../../output/05_qtl_analysis/01_liftover_vcf/tmp_files"
CONCATFILE="${SAVETMP}/tmp_concat_chroms.vcf"
TMPFILE="${SAVETMP}/tmp.vcf"
NORM_TMPFILE="${SAVETMP}/tmp_norm.vcf"

# Load necessary tools
echo "Loading required modules..."
module load gatk/4.0.10.0
module load bcftools

# Concatenate all of the chromosomes
echo "Concatenating files..."
bcftools concat "${ROSMAP_VCF_PATH}/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_{1..22}.recalibrated_variants.vcf.gz" \
    --output "${CONCATFILE}"

# Operations flags
STORE_DATA=true
NORMALIZE=true
LIFTOVER=true

# Add a header line for "chr" before chromosomes in the fasta file
if [ "$STORE_DATA" = true ]; then
    echo "Storing temp VCF..."
    cat \
        <(grep "^##" "$CONCATFILE") \
        <(cat "$HEADER_PATH") \
        <(grep -v "^##" "$CONCATFILE" | sed -r 's/^([0-9]+)/chr\1/') \
        > "${TMPFILE}"
fi

# Normalize the variant file to get more matches during liftover
if [ "$NORMALIZE" = true ]; then
    echo "Normalizing VCF..."
    bcftools norm \
        --fasta-ref="${FASTA_REF_HG19}" \
        --check-ref s \
        --multiallelics + \
        --output="${NORM_TMPFILE}" \
        "${TMPFILE}" 
fi

# Liftover from hg19 to hg38
if [ "$LIFTOVER" = true ]; then
    echo "Lifting over VCF..."
    gatk LiftoverVcf \
        --INPUT="${NORM_TMPFILE}" \
        --OUTPUT="${OUTPUT_FILE}" \
        --CHAIN="${CHAIN_FILE}" \
        --REFERENCE_SEQUENCE="${FASTA_REF_HG38}" \
        --REJECT="${REJECTED_FILE}" \
        --MAX_RECORDS_IN_RAM 5000 \
        --RECOVER_SWAPPED_REF_ALT
fi

# Final normalization
echo "Final normalization..."
bcftools norm \
    --fasta-ref="${FASTA_REF_HG38}" \
    --check-ref s \
    --multiallelics + \
    --output="${NORMALIZED_OUTPUT}" \
    "${OUTPUT_FILE}"

# End of Script

