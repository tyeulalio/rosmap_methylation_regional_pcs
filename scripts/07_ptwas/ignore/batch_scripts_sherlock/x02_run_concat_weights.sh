#!/usr/bin/env bash

#SBATCH --job-name=concat_ptwas_weights
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


SLURM_ARRAY_TASK_ID=0

# load variables from shared script
. ./helper_scripts/get_variables.sh

# select summary type here
SUMMARY_TYPE=${SUMMARY_TYPES[1]} # avgs, pcs
REGION_TYPE=${REGION_TYPES[2]} # full_gene, promoters, genome
CELL_TYPE=${CELL_TYPES[$SLURM_ARRAY_TASK_ID]} # bulk, astro, endo, neuron, oligo_opc

PROPORTION_TYPE="with_proportions"


# create the runname for this
. ./helper_scripts/get_runname.sh $SUMMARY_TYPE $CELL_TYPE $REGION_TYPE
echo $RUNNAME


OUTDIR="${HOME}/deconvolution/new_rosmap/output/07_ptwas/02_concatenated_weights"
mkdir -p $OUTDIR
OUTDIR="${OUTDIR}/${PROPORTION_TYPE}"
. ./helper_scripts/make_output_directory.sh $OUTDIR $RUNNAME

# directory containing the ptwas weights
WEIGHTSDIR="${HOME}/deconvolution/new_rosmap/output/07_ptwas/01_ptwas_weights/${PROPORTION_TYPE}/${RUNNAME}/ptwas_weights"


## --- run script here
../02_concat_weights.sh \
    $WEIGHTSDIR \
    $RUN_OUTDIR


# move the log files over
./helper_scripts/copy_logs.sh $SBATCH_PREFIX $LOG_OUTDIR $SLURM_ARRAY_TASK_ID
