#!/bin/bash

# using plink2 on scg
module load plink2

# run plink to generate files for rosmap


# make the plink2 format ped files
echo "Creating pgen files"
plink2 \
    --vcf "../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_all_chroms_normalized.vcf" \
    --make-pgen \
    --allow-extra-chr \
    --chr 1-22 \
    --output-chr chrM \
    --max-alleles 2 \
    --keep "../../output/01_preprocessing/01_subset_samples/keep_wgs_samples.tsv" \
    --mind 0.05 \
    --hwe 0.001 \
    --geno 0.05 \
    --maf 0.01 \
    --set-all-var-ids "@_#_\$r_\$a" \
    --new-id-max-allele-len 20 missing \
    --vcf-half-call missing \
    --out "../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms"


# make the plink format bed files
echo "Creating bed files"
plink2 \
    --vcf "../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_all_chroms_normalized.vcf" \
    --make-bed \
    --allow-extra-chr \
    --chr 1-22 \
    --output-chr chrM \
    --max-alleles 2 \
    --keep "../../output/01_preprocessing/01_subset_samples/keep_wgs_samples.tsv" \
    --mind 0.05 \
    --hwe 0.001 \
    --geno 0.05 \
    --maf 0.01 \
    --set-all-var-ids "@_#_\$r_\$a" \
    --new-id-max-allele-len 20 missing \
    --vcf-half-call missing \
    --out "../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms"


echo "Creating VCF files"
plink2 \
    --vcf "../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_all_chroms_normalized.vcf" \
    --keep "../../output/01_preprocessing/01_subset_samples/one_sample.tsv" \
    --export vcf \
    --allow-extra-chr \
    --chr 1-22 \
    --output-chr chrM \
    --max-alleles 2 \
    --set-all-var-ids "@_#_\$r_\$a" \
    --new-id-max-allele-len 20 missing \
    --vcf-half-call missing \
    --out "../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms"

    
    #--mind 0.05 \
    #--hwe 0.001 \
    #--geno 0.05 \
    #--maf 0.01 \
