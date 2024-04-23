#!/usr/bin/env bash

#SBATCH --job-name=ptwas_scan_allChroms
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --partition=nih_s10
#SBATCH --account=smontgom
#SBATCH --time=3-0:00:00
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
#REGION_TYPE=${REGION_TYPES[0]} # full_gene, promoters, genome
#CELL_TYPE=${CELL_TYPES[4]} # bulk, astro, endo, neuron, oligo_opc

SUMMARY_TYPE=$1
REGION_TYPE=$2
CELL_TYPE=$3

PROPORTION_TYPE="with_proportions"
GWAS_TYPE="wightman"

# create the runname for this
. ./helper_scripts/get_runname.sh $SUMMARY_TYPE $CELL_TYPE $REGION_TYPE
echo $RUNNAME


OUTDIR="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/07_ptwas_scan/${GWAS_TYPE}"
mkdir -p $OUTDIR
OUTDIR="${OUTDIR}/${PROPORTION_TYPE}"
. ./helper_scripts/make_output_directory.sh $OUTDIR $RUNNAME


# GWAS TYPE
GWAS_FILE="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/04_formatted_gwas/${GWAS_TYPE}_gwas_file.vcf.gz"

#PTWAS_WEIGHT_FILE="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/03_gambit_db/${RUNNAME}/formatted_all_gene.ptwas_weights.gambit.vcf.gz"
CHROM=$SLURM_ARRAY_TASK_ID
#CHROM=22

PTWAS_WEIGHT_FILE="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/03_gambit_db/${PROPORTION_TYPE}/${RUNNAME}/all_gene.ptwas_weights.gambit.vcf.gz"

NUM_LINES=$(zcat ${PTWAS_WEIGHT_FILE} | wc -l)
echo "num lines ${NUM_LINES}"

if [ $NUM_LINES -lt 3 ]; then
    echo "Less than 3 lines in input file. Skipping"
    # move the log files over
    ./helper_scripts/copy_logs.sh $SBATCH_PREFIX $LOG_OUTDIR $SLURM_ARRAY_TASK_ID
    exit
fi

LD_PANEL_FILES="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/06_swapped_ld_panels/$GWAS_TYPE/chr*.vcf.gz"
set -f


OUTPUT_PREFIX="${RUN_OUTDIR}/all_gene_ptwas_scan"

## --- run script here
../07_ptwas_scan.sh \
    $GWAS_FILE \
    $PTWAS_WEIGHT_FILE \
    $LD_PANEL_FILES \
    $OUTPUT_PREFIX


# move the log files over
./helper_scripts/copy_logs.sh $SBATCH_PREFIX $LOG_OUTDIR $SLURM_ARRAY_TASK_ID
