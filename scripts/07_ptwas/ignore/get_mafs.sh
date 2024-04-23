#!/usr/bin/env bash

# get the MAFs from ROSMAP VCFs formatted nicely


# make sure conda environment is activated
. ./batch_scripts/helper_scripts/activate_colocalization_env.sh


# creat the savedir
SAVEDIR="/home/eulalio/deconvolution/new_rosmap/data/rosmap_wgs_liftover_hg38/MAFs"
echo "Writing output to $SAVEDIR"
mkdir -p $SAVEDIR

#CHROMFILE="/home/eulalio/deconvolution/new_rosmap/data/rosmap_wgs_liftover_hg38/rosmap_wgs_hg38_chr${CHROM}.vcf"
VCFFILE="/home/eulalio/deconvolution/new_rosmap/data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms.vcf.gz"
OUTFILE="${SAVEDIR}/rosmap_wgs_hg38_all_chroms_MAFs.txt"

bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" $VCFFILE > $OUTFILE

