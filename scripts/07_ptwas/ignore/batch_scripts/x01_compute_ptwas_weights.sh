#!/usr/bin/env bash

#SBATCH --job-name=ptwas_weights
#SBATCH --array=0-499
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=nih_s10
#SBATCH --account=smontgom
#SBATCH --time=0-12:00:00
#SBATCH --wait
#SBATCH --output=/home/eulalio/deconvolution/new_rosmap/output/sbatch_output/%x_array_%a.out
#SBATCH --error=/home/eulalio/deconvolution/new_rosmap/output/sbatch_output/%x_array_%a.err
SBATCH_PREFIX="/home/eulalio/deconvolution/new_rosmap/output/sbatch_output/${SLURM_JOB_NAME}_array"

# make sure conda environment is activated
. ./helper_scripts/activate_colocalization_env.sh

# load variables from shared script
. ./helper_scripts/get_variables.sh

# select summary type here
#SUMMARY_TYPE=${SUMMARY_TYPES[0]} # avgs, pcs
#REGION_TYPE=${REGION_TYPES[8]} # full_gene, promoters, genome
#CELL_TYPE=${CELL_TYPES[0]} # bulk, astro, endo, neuron, oligo_opc
SUMMARY_TYPE=$1
REGION_TYPE=$2
CELL_TYPE=$3
TOTAL_GENES=$4

#SLURM_ARRAY_TASK_ID=1

PROPORTION_TYPE="with_proportions"

TOTAL_JOBS=499


echo "Processing ${CELL_TYPE}"

# create the runname for this
. ./helper_scripts/get_runname.sh $SUMMARY_TYPE $CELL_TYPE $REGION_TYPE
echo $RUNNAME


OUTDIR="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/01_ptwas_weights"
mkdir -p $OUTDIR
OUTDIR="${OUTDIR}/${PROPORTION_TYPE}"
. ./helper_scripts/make_output_directory.sh $OUTDIR $RUNNAME

# input directories
SBAMSDIR="/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/03_dapg_input/${PROPORTION_TYPE}/${RUNNAME}/sbams"
DAPDIR="/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/04_dapg/${PROPORTION_TYPE}/${RUNNAME}"


# get list of genes 
GENES_LIST=$(ls $DAPDIR | sed -r 's/.dapg//g')
FILE_COUNT=$(echo $GENES_LIST | wc -w)
echo "Total files: $FILE_COUNT"

# create array of files
readarray -t GENES_ARRAY <<<"$GENES_LIST"

# find start and end for this index
result=$(echo "scale=2; $TOTAL_GENES / $TOTAL_JOBS" | bc)
ceiling=$(echo "scale=0; ($result+1)/1" | bc)
echo "Result: $result"
echo "Ceiling: $ceiling"

GENES_PER_JOB=$ceiling

START_IDX=$(((SLURM_ARRAY_TASK_ID-1) * GENES_PER_JOB + 1))
END_IDX=$((START_IDX + GENES_PER_JOB))

echo "id $SLURM_ARRAY_TASK_ID"
echo "start $START_IDX"
echo "end $END_IDX"

MIN_IDX=0

if [ "$START_IDX" -lt "$MIN_IDX" ]; then
    START_IDX=$MIN_IDX
fi

if [ "$END_IDX" -gt "$TOTAL_GENES" ]; then
    END_IDX=$TOTAL_GENES
fi

# make sure we do gene 0
if [ "$START_IDX" -eq 1 ]; then
    START_IDX=0
fi

echo "start $START_IDX"
echo "end $END_IDX"


time for i in $(seq $START_IDX $END_IDX); do
    echo "starting $i"
    GENE_NUM=$i
    #GENE_NUM=1

    # select the gene that we're working with
    GENE=${GENES_ARRAY[$GENE_NUM]}

    # run ptwas builder for each gene
    # file containing sbam info
    SBAMGENE=$(echo $GENE | sed 's/-/./1')
    SBAMFILE="${SBAMSDIR}/${SBAMGENE}_fastqtl_singl_snp_output.sbam"
    echo "SBAM FILE: ${SBAMFILE}"

    # file containing dap-g output
    DAPGENE=$(echo $GENE | tr \\. -)
    DAPFILE="${DAPDIR}/${DAPGENE}.dapg"
    echo "DAP FILE: ${DAPFILE}"


    echo "processing $GENE $GENE_NUM out of $FILE_COUNT"
    ((i+=1))

    ## --- run script here
    ../01_compute_ptwas_weights.sh \
        $DAPFILE \
        $SBAMFILE \
        $GENE \
        $RUN_OUTDIR
    wait
done


# move the log files over
./helper_scripts/copy_logs.sh $SBATCH_PREFIX $LOG_OUTDIR $SLURM_ARRAY_TASK_ID
