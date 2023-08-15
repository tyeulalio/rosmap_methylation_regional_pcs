#!/bin/bash

# liftover the VCFs from hg17 to hg19
# using GATK liftoverVCF


# make sure to load the GATK module
module load gatk/4.0.10.0
module load bcftools

#CHROM=$1
#CHROM=22


#for CHROM in {22..8}
#do
    SAVETMP="../../output/05_qtl_analysis/01_liftover_vcf/tmp_files"
    #SAVETMP=${TMPDIR}

    # storing temporary files for the run
    CONCATFILE="${SAVETMP}/tmp_concat_chroms.vcf"
    TMPFILE="${SAVETMP}/tmp.vcf"
    NORM_TMPFILE="${SAVETMP}/tmp_norm.vcf"
    #FORM_TMPFILE="${SAVETMP}/tmp_formatted.vcf"

    # concatenate all of the chromosomes first
    #echo "Concatenating files"
    #bcftools concat /home/eulalio/deconvolution/rosmap/data/rosmap_wgs_harmonization_variant_calling/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_{1..22}.recalibrated_variants.vcf.gz \
        #--output ${CONCATFILE}


    # select which functions to perform
    STORE_DATA=false
    NORMALIZE=false
    LIFTOVER=false

    # move to tmp folder and unzip
    # adds a header line to accounts for "chr" before chromosomes in the fasta file
    # adds "chr" before the chromosome in the wgs vcf
    if [ "$STORE_DATA" = true ]
    then
        echo "Storing temp VCF"
        cat \
            <(cat "$CONCATFILE" | grep "^##" ) \
            <(cat "/home/eulalio/deconvolution/data/reference/vcf_header_chr_lines.txt") \
            <(cat "$CONCATFILE" | grep -v "^##" | sed -r 's/^([0-9]+)/chr\1/') \
            > ${TMPFILE}
    fi

    # normalize the variant file first
    # can help to get more matches during liftOver
    if [ "$NORMALIZE" = true ]
    then
        echo "Normalizing vcf"
        bcftools norm \
            --fasta-ref="/reference/RefGenomes/GATK_Resource_Bundle/hg19/hg19.fa" \
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
        --OUTPUT="../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_hg38_all_chroms.vcf" \
        --CHAIN="/home/eulalio/shared/liftOver/chains/hg19ToHg38.over.chain.gz" \
        --REFERENCE_SEQUENCE="/reference/RefGenomes/GATK_Resource_Bundle/hg38/hg38.fa" \
        --REJECT="../../output/05_qtl_analysis/01_liftover_vcf/rejected_records_all_chroms.vcf" \
        --MAX_RECORDS_IN_RAM 5000 \
        --RECOVER_SWAPPED_REF_ALT
        #--VERBOSITY=DEBUG
    fi


        echo "Normalizing vcf"
        bcftools norm \
            --fasta-ref="/reference/RefGenomes/GATK_Resource_Bundle/hg38/hg38.fa" \
            --check-ref s \
            --multiallelics + \
            --output="../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_all_chroms_normalized.vcf" \
            "../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_hg38_all_chroms.vcf" 

    exit


    echo "Processing chromosome ${CHROM}"

    # select which functions to perform
    STORE_DATA=true
    NORMALIZE=true
    LIFTOVER=true


    # move to tmp folder and unzip
    # adds a header line to accounts for "chr" before chromosomes in the fasta file
    # adds "chr" before the chromosome in the wgs vcf
    if [ "$STORE_DATA" = true ]
    then
        echo "Storing temp VCF"
        cat \
            <(zcat "/home/eulalio/deconvolution/rosmap/data/rosmap_wgs_harmonization_variant_calling/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_${CHROM}.recalibrated_variants.vcf.gz" | grep "^##" ) \
            <(cat "/home/eulalio/deconvolution/data/reference/vcf_header_chr_lines.txt") \
            <(zcat "/home/eulalio/deconvolution/rosmap/data/rosmap_wgs_harmonization_variant_calling/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_${CHROM}.recalibrated_variants.vcf.gz" | grep -v "^##" | sed -r 's/^([0-9]+)/chr\1/') \
            > ${TMPFILE}
    fi


    # normalize the variant file first
    # can help to get more matches during liftOver
    if [ "$NORMALIZE" = true ]
    then
        echo "Normalizing vcf"
        bcftools norm \
            --fasta-ref="/reference/RefGenomes/GATK_Resource_Bundle/hg19/hg19.fa" \
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
        --OUTPUT="../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_hg38_chr${CHROM}.vcf" \
        --CHAIN="/home/eulalio/shared/liftOver/chains/hg19ToHg38.over.chain.gz" \
        --REFERENCE_SEQUENCE="/reference/RefGenomes/GATK_Resource_Bundle/hg38/hg38.fa" \
        --REJECT="../../output/05_qtl_analysis/01_liftover_vcf/rejected_records_chr${CHROM}.vcf" \
        --MAX_RECORDS_IN_RAM 5000 \
        --RECOVER_SWAPPED_REF_ALT
        #--VERBOSITY=DEBUG
    fi

    # remove the temp file
    #rm "$TMPFILE"

#done

        #--REFERENCE_SEQUENCE="/home/eulalio/shared/genomes/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa" \
        #--INPUT="/home/eulalio/deconvolution/rosmap/data/rosmap_wgs_harmonization_variant_calling/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_${CHROM}.recalibrated_variants.vcf.gz" \
