#!/bin/bash

# Scripts to liftover genotype VCFs from hg17 to hg19 using GATK liftoverVCF

# load required modules
module load gatk/4.0.10.0
module load bcftools

# Parameters
# ROSMAP_VCF_PATH = directory containing rosmap VCF files
# HEADER_PATH = VCF header to use at the top of each file, just take from a valid VCF
# FASTA_REF_HG19 = fasta reference file for hg19 build
# OUTPUT_FILE = path to save the final lifted output
# CHAIN_FILE = liftOver chain file, can be downloaded from UCSC
# FASTA_REF_HG19 = fasta reference file for hg38 build
# REJECTED_FILE = file to write SNPs that failed liftOver
# NORMALIZED_OUTPUT_FILE = normalized, lifted output will be saved here
ROSMAP_VCF_PATH="/path/to/rosmap_wgs_harmonization_variant_calling"
HEADER_PATH="/home/eulalio/deconvolution/data/reference/vcf_header_chr_lines.txt"
FASTA_REF_HG19="/reference/RefGenomes/GATK_Resource_Bundle/hg19/hg19.fa" 
OUTPUT_FILE="../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_hg38_all_chroms.vcf" 
CHAIN_FILE="/home/eulalio/shared/liftOver/chains/hg19ToHg38.over.chain.gz" 
FASTA_REF_HG38="/reference/RefGenomes/GATK_Resource_Bundle/hg38/hg38.fa" 
REJECTED_FILE="../../output/05_qtl_analysis/01_liftover_vcf/rejected_records_all_chroms.vcf" 
NORMALIZED_OUTPUT_FILE="../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_all_chroms_normalized.vcf" 

# Temporary storage for operations
SAVETMP="../../output/05_qtl_analysis/01_liftover_vcf/tmp_files"
CONCATFILE="${SAVETMP}/tmp_concat_chroms.vcf"
TMPFILE="${SAVETMP}/tmp.vcf"
NORM_TMPFILE="${SAVETMP}/tmp_norm.vcf"



# concatenate all of the chromosomes first
echo "Concatenating files"
bcftools concat ${ROSMAP_VCF_PATH}/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_{1..22}.recalibrated_variants.vcf.gz \
    --output ${CONCATFILE}


# select which functions to perform
STORE_DATA=true
NORMALIZE=true
LIFTOVER=true

# adds a header line to accounts for "chr" before chromosomes in the fasta file
# adds "chr" before the chromosome in the wgs vcf
if [ "$STORE_DATA" = true ]
then
    echo "Storing temp VCF"
    cat \
        <(cat "$CONCATFILE" | grep "^##" ) \
        <(cat "$HEADERPATH") \
        <(cat "$CONCATFILE" | grep -v "^##" | sed -r 's/^([0-9]+)/chr\1/') \
        > ${TMPFILE}
fi


# normalize the variant file first
# can help to get more matches during liftOver
if [ "$NORMALIZE" = true ]
then
    echo "Normalizing vcf"
    bcftools norm \
        --fasta-ref=${FASTA_REF_HG19} \
        --check-ref s \
        --multiallelics + \
        --output="$NORM_TMPFILE" \
        "${TMPFILE}" 
fi


# perform liftover
# liftover from hg19 to hg38
if [ "$LIFTOVER" = true ]
then
echo "Lifting over vcf"
gatk LiftoverVcf \
    --INPUT="$NORM_TMPFILE" \
    --OUTPUT="$OUTPUT_FILE" \
    --CHAIN="$CHAIN_FILE" \
    --REFERENCE_SEQUENCE="$FASTA_REF_HG38" \
    --REJECT="$REJECTED_FILE" \
    --MAX_RECORDS_IN_RAM 5000 \
    --RECOVER_SWAPPED_REF_ALT
fi


# final normalization
echo "Final normalization"
bcftools norm \
    --fasta-ref="$FASTA_REF_HG38" \
    --check-ref s \
    --multiallelics + \
    --output="$NORMALIZED_OUTPUT" \
    "$OUTPUT_FILE"


