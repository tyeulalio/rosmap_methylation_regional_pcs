#!/home/eulalio/micromamba/bin/bash


# split the normalized rosmap wgs vcf

CHROM=$1

# rosmap input file
#ROSMAP_FILE="/home/eulalio/deconvolution/new_rosmap/data/rosmap_wgs_liftover_hg38/rosmap_wgs_all_chroms_normalized.vcf"
INFILE="/home/eulalio/deconvolution/new_rosmap/data/rosmap_wgs_liftover_hg38/rosmap_wgs_hg38_chr${CHROM}.vcf"

# output file for chromosome
OUTDIR="/home/eulalio/deconvolution/new_rosmap/data/rosmap_wgs_liftover_hg38/normalized_vcfs"
mkdir -p $OUTDIR

OUTFILE="${OUTDIR}/rosmap_hgs_hg38_normalized_chrom${CHROM}.vcf"

# try bcftools
bcftools \
    annotate -x INFO,^FORMAT/GT \
    $INFILE \
    -O v \
    -o /tmp/tmp.txt

    #-r "$CHROM" \
#plink \
    #--vcf $ROSMAP_FILE \
    #--chr $CHROM \
    #--export vcf \
    #--out $OUTFILE
    
