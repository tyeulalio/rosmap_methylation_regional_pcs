#!/usr/bin/env bash

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate colocalization

# run plink to generate files for rosmap
# need these gwas vcfs to create LD matrix for DAP
# this only needs to be run once - do not repeat for every summary type/region type/ cell type

OUTDIR="../../../output/06_fine_mapping_colocalization/02_gwas_prep/03_create_gwas_vcfs/"

mkdir -p ${OUTDIR}

PLINKDIR="${OUTDIR}plink/"
LDDIR="${OUTDIR}ld/"

# create vcf with all chromosomes
# takes a long time, just run once and comment out
plink2 \
    --vcf "../../data/rosmap_wgs_liftover_hg38/rosmap_wgs_all_chroms_normalized.vcf" \
    --make-bed \
    --allow-extra-chr \
    --chr 1-22 \
    --output-chr chrM \
    --mind 0.05 \
    --geno 0.05 \
    --maf 0.01 \
    --max-alleles 2 \
    --set-all-var-ids "@_#_\$r_\$a" \
    --new-id-max-allele-len 20 missing \
    --vcf-half-call missing \
    --out "${PLINKDIR}rosmap_wgs_hg38_all_chroms" \
    > >(tee -a ${PLINKDIR}plink2.stdout 2> >(tee -a ${PLINKDIR}plink2.stderr >&2))

# check on snps with ID .
plink2 \
    --export vcf \
    --snp "." \
    --bfile "${PLINKDIR}rosmap_wgs_hg38_all_chroms" \
    --out "${PLINKDIR}rosmap_wgs_hg38_all_chroms_dup_snp" 

# remove duplicate IDs
plink2 \
    --make-bed \
    --exclude-snp "." \
    --bfile "${PLINKDIR}rosmap_wgs_hg38_all_chroms" \
    --out "${PLINKDIR}rosmap_wgs_hg38_all_chroms_no_dups" 
