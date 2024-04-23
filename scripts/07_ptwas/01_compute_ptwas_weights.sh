#!/usr/bin/env bash

# command line input
DAPFILE=$1
SBAMFILE=$2
GENE=$3
OUTDIR=$4

# format output file
# put genes in same format
GENE=$(echo $GENE | tr - .)

WEIGHTS_OUTDIR="${OUTDIR}/ptwas_weights/"
mkdir -p ${WEIGHTS_OUTDIR}
OUTPUTFILE="${WEIGHTS_OUTDIR}/${GENE}_ptwas_weights.txt"

# run ptwas bulider
ptwas_builder \
    -f $DAPFILE \
    -d $SBAMFILE \
    -g $GENE \
    -o $OUTPUTFILE
