#!/bin/bash

# using plink2 on scg
module load plink2

# run plink to generate files for rosmap


# -----------------
# PARAMETER OPTIONS
# -----------------

INPUT_VCF_PATH="../path_to_your_data/input_vcf_file.vcf"
OUTPUT_PREFIX="../path_to_output_dir/output_prefix"

SAMPLES_KEEP_PATH="../path_to_samples_dir/samples_to_keep.tsv"
ONE_SAMPLE_PATH="../path_to_samples_dir/one_sample.tsv"

# Plink2 Filters (adjust as needed)
MIND_FILTER=0.05
HWE_FILTER=0.001
GENO_FILTER=0.05
MAF_FILTER=0.01

# Load necessary modules (if using environment modules system)
module load plink2

# Convert VCF to PLINK's PGEN format
echo "Creating pgen files"
plink2 \
    --vcf "$INPUT_VCF_PATH" \
    --make-pgen \
    --allow-extra-chr \
    --chr 1-22 \
    --output-chr chrM \
    --max-alleles 2 \
    --keep "$SAMPLES_KEEP_PATH" \
    --mind $MIND_FILTER \
    --hwe $HWE_FILTER \
    --geno $GENO_FILTER \
    --maf $MAF_FILTER \
    --set-all-var-ids "@_#_\$r_\$a" \
    --new-id-max-allele-len 20 missing \
    --vcf-half-call missing \
    --out "${OUTPUT_PREFIX}_pgen"


# Convert VCF to PLINK's BED format
echo "Creating bed files"
plink2 \
    --vcf "$INPUT_VCF_PATH" \
    --make-bed \
    --allow-extra-chr \
    --chr 1-22 \
    --output-chr chrM \
    --max-alleles 2 \
    --keep "$SAMPLES_KEEP_PATH" \
    --mind $MIND_FILTER \
    --hwe $HWE_FILTER \
    --geno $GENO_FILTER \
    --maf $MAF_FILTER \
    --set-all-var-ids "@_#_\$r_\$a" \
    --new-id-max-allele-len 20 missing \
    --vcf-half-call missing \
    --out "${OUTPUT_PREFIX}_bed"

# Extract a sample from the VCF
echo "Creating VCF files for a single sample"
plink2 \
    --vcf "$INPUT_VCF_PATH" \
    --keep "$ONE_SAMPLE_PATH" \
    --export vcf \
    --allow-extra-chr \
    --chr 1-22 \
    --output-chr chrM \
    --max-alleles 2 \
    --set-all-var-ids "@_#_\$r_\$a" \
    --new-id-max-allele-len 20 missing \
    --vcf-half-call missing \
    --out "${OUTPUT_PREFIX}_one_sample_vcf"
