#!/usr/bin/env bash

# call the ptwas_est step

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate colocalization

PTWAS_INPUT_FILE=$1
GENE_NAME=$2
SAVEFILE=$3
PIP_THRESH=$4

echo "saving to $SAVEFILE"

PTWAS_est -d $PTWAS_INPUT_FILE \
    -t $PIP_THRESH \
    -n $GENE_NAME > \
    "$SAVEFILE"
