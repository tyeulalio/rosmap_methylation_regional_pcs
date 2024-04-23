#!/usr/bin/env bash

#SBATCH --job-name=ptwas_scan_sher
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --account=smontgom
#SBATCH --time=1-00:00:00
#SBATCH --wait
#SBATCH --output=/home/users/eulalio/deconvolution/new_rosmap/output/sbatch_output/%x_array_%a.out
#SBATCH --error=/home/users/eulalio/deconvolution/new_rosmap/output/sbatch_output/%x_array_%a.err
SBATCH_PREFIX="/home/users/eulalio/deconvolution/new_rosmap/output/sbatch_output/${SLURM_JOB_NAME}_array"

# make sure conda environment is activated
# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate /oak/stanford/groups/smontgom/eulalio/micromamba/colocalization

# load variables from shared script
. ./helper_scripts/get_variables.sh

# select summary type here
#SUMMARY_TYPE=${SUMMARY_TYPES[1]} # avgs, pcs
#REGION_TYPE=${REGION_TYPES[0]} # full_gene, promoters, genome
#CELL_TYPE=${CELL_TYPES[2]} # bulk, astro, endo, neuron, oligo_opc
PROPORTION_TYPE="with_proportions"

SUMMARY_TYPE=$1
REGION_TYPE=$2
CELL_TYPE=$3

#SLURM_ARRAY_TASK_ID=22
CHROM=$SLURM_ARRAY_TASK_ID

GWAS_TYPE="wightman"

# create the runname for this
. ./helper_scripts/get_runname.sh $SUMMARY_TYPE $CELL_TYPE $REGION_TYPE
echo $RUNNAME


OUTDIR="${HOME}/deconvolution/new_rosmap/output/07_ptwas/07_ptwas_scan/${GWAS_TYPE}"
mkdir -p $OUTDIR
OUTDIR="${OUTDIR}/${PROPORTION_TYPE}"
. ./helper_scripts/make_output_directory.sh $OUTDIR $RUNNAME


# GWAS TYPE
GWAS_FILE="${HOME}/deconvolution/new_rosmap/output/07_ptwas/04_formatted_gwas/${GWAS_TYPE}_gwas_file.vcf.gz"

#PTWAS_WEIGHT_FILE="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/03_gambit_db/${RUNNAME}/formatted_all_gene.ptwas_weights.gambit.vcf.gz"

PTWAS_WEIGHT_FILE="${HOME}/deconvolution/new_rosmap/output/07_ptwas/03_gambit_db/${PROPORTION_TYPE}/${RUNNAME}/chrom${CHROM}_all_gene.ptwas_weights.gambit.vcf.gz"

NUM_LINES=$(zcat ${PTWAS_WEIGHT_FILE} | wc -l)
echo "num lines ${NUM_LINES}"

if [ $NUM_LINES -lt 3 ]; then
    echo "Less than 3 lines in input file. Skipping"
    # move the log files over
    ./helper_scripts/copy_logs.sh $SBATCH_PREFIX $LOG_OUTDIR $SLURM_ARRAY_TASK_ID
    exit
fi

LD_PANEL_FILES="${HOME}/deconvolution/new_rosmap/output/07_ptwas/06_swapped_ld_panels/$GWAS_TYPE/chr*.vcf.gz"
set -f


OUTPUT_PREFIX="${RUN_OUTDIR}/chr${CHROM}_ptwas_scan"

## --- run script here
../07_ptwas_scan.sh \
    $GWAS_FILE \
    $PTWAS_WEIGHT_FILE \
    $LD_PANEL_FILES \
    $OUTPUT_PREFIX


# move the log files over
./helper_scripts/copy_logs.sh $SBATCH_PREFIX $LOG_OUTDIR $SLURM_ARRAY_TASK_ID
