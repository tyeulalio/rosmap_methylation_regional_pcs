#!/usr/bin/env bash

# get a template file for the genotype counts for each genotype
# the templates are used by script 03


# look into sbam files that have already been generated by 03 script
SBAM_DIR="/path/to/output/06_fine_mapping_colocalization/03_qtl_prep/03_dapg_input/without_proportions/avgs_astro_full_gene/sbams"

SBAM_FILES=$(ls $SBAM_DIR)

# output directory
OUTPUT_DIR="/path/to/output/06_fine_mapping_colocalization/03_qtl_prep/sbam_genotype_files"
mkdir -p $OUTPUT_DIR

i=0
for SBAM_FILE in $SBAM_FILES
do
    #echo $SBAM_FILE

    # get gene from sbam file name
    GENE=${SBAM_FILE%%_*}
    #echo $GENE

    OUTFILE="${OUTPUT_DIR}/${GENE}_genotype_counts.txt"
    #echo "write output to $OUTFILE"

    grep 'geno ' $SBAM_DIR/$SBAM_FILE > $OUTFILE

    ((i++))
    if [ $((i % 1000)) -eq 0 ];then
        echo "Processed $i genes"
    fi

done