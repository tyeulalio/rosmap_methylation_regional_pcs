#!/usr/bin/env bash

# load conda environment
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
#mamba activate base
mamba activate r

PTWAS_WEIGHTS_FILE=$1
OUTDIR=$2

# output file
OUTPUT_FILE="${OUTDIR}/all_gene.ptwas_weights.gambit.vcf"
rm -f $OUTPUT_FILE
touch $OUTPUT_FILE

# get command line input
# run the r script
Rscript ../make_GAMBIT.DB.R \
    -d ${PTWAS_WEIGHTS_FILE} \
    -o ${OUTPUT_FILE}

# format the output file
mamba activate colocalization

# format the output file by chromsome
CHROM=1


for CHROM in {1..22};
do
    echo "Processing chromosome ${CHROM}"

    # get the header from gabmit file
    FORMATTED_OUTFILE="${OUTDIR}/chrom${CHROM}_all_gene.ptwas_weights.gambit.vcf"
    touch $FORMATTED_OUTFILE

    echo "Writing output to ${FORMATTED_OUTFILE}"

    echo "printing header"
    cat $OUTPUT_FILE | head | grep "^#" > $FORMATTED_OUTFILE

    # remove chr from chromosome and from ID
    echo "printing remaining formatted lines"
    cat $OUTPUT_FILE | grep -v "^#" | grep "^chr${CHROM}\s" | sed -e "s/chr//g" >> $FORMATTED_OUTFILE

    # bgzip the output file
    echo "zipping output file"
    bgzip -f $FORMATTED_OUTFILE

    # tabix the output
    echo "running tabix"
    tabix -f -p 'vcf' "${FORMATTED_OUTFILE}.gz"
done
