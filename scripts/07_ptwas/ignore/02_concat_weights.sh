#!/usr/bin/env bash

# concatenate the ptwas weights 

WEIGHTSDIR=$1
OUTDIR=$2

PTWAS_FILES=$(ls ${WEIGHTSDIR} | grep '_ptwas_weights')

# creat the outputfile
OUTFILE="${OUTDIR}/all_gene.ptwas_weights.txt"
rm -f $OUTFILE
touch $OUTFILE

for PTWAS_FILE in $PTWAS_FILES; do
    cat "${WEIGHTSDIR}/${PTWAS_FILE}" >> $OUTFILE
done
#cat ${WEIGHTSDIR}/*_ptwas_weights.txt | gzip - > "${OUTDIR}/all_gene.ptwas_weights.gz"
#

# remove final file if it exists
ZIP_OUTFILE="${OUTDIR}/all_gene.ptwas_weights.txt.gz"
rm -f $ZIP_OUTFILE

# gzip the file
gzip -f $OUTFILE 
