#!/usr/bin/env bash

# check counts for results

source ./batch_scripts/helper_scripts/get_variables.sh


SUMMARY_TYPE='avgs'
#SUMMARY_TYPE='pcs'
REGION_TYPE='preTSS'
#CELL_TYPE='oligo_opc'

CHECK01=false
CHECK02=false
CHECK03=true

# check 01 results
if [ "$CHECK01" = true ]; then
    for CELL_TYPE in "${CELL_TYPES[@]}"
    do
        echo $CELL_TYPE

        echo "Checking 01"
        RUNNAME="${SUMMARY_TYPE}_${CELL_TYPE}_${REGION_TYPE}"
        DIR01="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/01_ptwas_weights/with_proportions/${RUNNAME}/ptwas_weights"
        echo $DIR01
        ls $DIR01 | wc -l


    done
fi

# check 02 script
if [ "$CHECK02" = true ]; then
    for CELL_TYPE in "${CELL_TYPES[@]}"
    do
        echo $CELL_TYPE

        echo "Checking 02"
        RUNNAME="${SUMMARY_TYPE}_${CELL_TYPE}_${REGION_TYPE}"

        DIR02="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/02_concatenated_weights/with_proportions/${RUNNAME}"
        echo $DIR02
        zcat "$DIR02/all_gene.ptwas_weights.txt.gz" | cut -f2 | uniq | wc -l

    done
fi

# check 03 script
if [ "$CHECK03" = true ]; then
    for CELL_TYPE in "${CELL_TYPES[@]}"
    do
        echo $CELL_TYPE

        echo "Checking 03"
        RUNNAME="${SUMMARY_TYPE}_${CELL_TYPE}_${REGION_TYPE}"

        DIR03="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/03_gambit_db/with_proportions/${RUNNAME}"
        echo $DIR03

        for CHROM in 1 5 10 15 22
        do
            echo "$RUNNAME -- Chromosome $CHROM"
            zcat "$DIR03/chrom${CHROM}_all_gene.ptwas_weights.gambit.vcf.gz" | wc -l
        done

    done
fi
