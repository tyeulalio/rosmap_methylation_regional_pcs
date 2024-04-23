#!/usr/bin/env bash

# run x01 script for each cell type
. ./helper_scripts/get_variables.sh

SUMMARY_TYPE_NUM=1
REGION_TYPE_NUM=3



SUMMARY_TYPE=${SUMMARY_TYPES[$SUMMARY_TYPE_NUM]} # avgs, pcs
REGION_TYPE=${REGION_TYPES[$REGION_TYPE_NUM]} # full_gene, promoters, exons

# run for each cell type
for CELL_TYPE_NUM in {0..4}
do
    CELL_TYPE=${CELL_TYPES[$CELL_TYPE_NUM]} # bulk, astro, endo, neuron, oligo_opc

    echo "$SUMMARY_TYPE $REGION_TYPE $CELL_TYPE"

    sbatch ./x07_run_ptwas_scan.sh $SUMMARY_TYPE $REGION_TYPE $CELL_TYPE 
    wait
done
