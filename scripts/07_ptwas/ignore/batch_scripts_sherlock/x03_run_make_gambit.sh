#!/usr/bin/env bash

#SBATCH --job-name=make_gambit
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --account=smontgom
#SBATCH --time=0-6:00:00
#SBATCH --wait
#SBATCH --output=/home/users/eulalio/deconvolution/new_rosmap/output/sbatch_output/%x_array_%a.out
#SBATCH --error=/home/users/eulalio/deconvolution/new_rosmap/output/sbatch_output/%x_array_%a.err
SBATCH_PREFIX="${HOME}/deconvolution/new_rosmap/output/sbatch_output/${SLURM_JOB_NAME}_array"


# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate /oak/stanford/groups/smontgom/eulalio/micromamba/colocalization

# load variables from shared script
. ./helper_scripts/get_variables.sh

SLURM_ARRAY_TASK_ID=1

# select summary type here
SUMMARY_TYPE=${SUMMARY_TYPES[1]} # avgs, pcs
REGION_TYPE=${REGION_TYPES[0]} # full_gene, promoters, genome
CELL_TYPE=${CELL_TYPES[$SLURM_ARRAY_TASK_ID]} # bulk, astro, endo, neuron, oligo_opc
#CELL_TYPE=${CELL_TYPES[1]} # bulk, astro, endo, neuron, oligo_opc
PROPORTION_TYPE="with_proportions"


# create the runname for this
. ./helper_scripts/get_runname.sh $SUMMARY_TYPE $CELL_TYPE $REGION_TYPE
echo $RUNNAME


OUTDIR="${HOME}/deconvolution/new_rosmap/output/07_ptwas/03_gambit_db"
mkdir -p $OUTDIR
OUTDIR="${OUTDIR}/${PROPORTION_TYPE}"
. ./helper_scripts/make_output_directory.sh $OUTDIR $RUNNAME


PTWAS_WEIGHTS_FILE="${HOME}/deconvolution/new_rosmap/output/07_ptwas/02_concatenated_weights/${PROPORTION_TYPE}/${RUNNAME}/all_gene.ptwas_weights.txt.gz"


## --- run script here
../03_run_make_gambit.sh \
    $PTWAS_WEIGHTS_FILE \
    $RUN_OUTDIR


# move the log files over
./helper_scripts/copy_logs.sh $SBATCH_PREFIX $LOG_OUTDIR $SLURM_ARRAY_TASK_ID
