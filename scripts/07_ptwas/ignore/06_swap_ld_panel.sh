#!/home/eulalio/micromamba/bin/bash

# swap reference/alt alleles for subset of variants in the LD panel
# to match the GWAS and QTL annotations

#CHROM=$1
GWAS_TYPE="wightman"

for CHROM in {22..21}
do
    echo "Processing chromosome ${CHROM}"
    # use the lifted LD panel
    LD_VCF="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/05_lifted_ld_panel/${GWAS_TYPE}/G1K_EUR_3V5_hg38/chr${CHROM}.vcf.gz"

    # swapped file containing alelles to update
    SWAPFILE="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/05_lifted_ld_panel/${GWAS_TYPE}/plink_formatted_swapped_chr${CHROM}.txt"

    # output file
    OUTDIR="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/06_swapped_ld_panels/${GWAS_TYPE}"
    mkdir -p $OUTDIR
    OUTFILE="${OUTDIR}/chr${CHROM}"

    TMPFILE="/tmp/tmp.txt"

    # run plink
    plink \
        --vcf $LD_VCF \
        --update-alleles $SWAPFILE \
        --export vcf \
        --out $OUTFILE

    # change the output
    echo "formatting genotype data"
    cat ${OUTFILE}.vcf | sed 's/\([0-9]\)\/\([0-9]\)/\1|\2/g' > $TMPFILE

    # bgzip the output file
    echo "bgzipping file"
    bgzip -f -c "${TMPFILE}" > "${OUTFILE}.vcf.gz"

    # try running tabix
    echo "creating file index"
    tabix -f -p vcf "${OUTFILE}.vcf.gz"
done
