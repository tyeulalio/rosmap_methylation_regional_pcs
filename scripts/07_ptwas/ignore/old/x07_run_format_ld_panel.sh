#!/usr/bin/env bash

#SBATCH --job-name=format_ld_panel
#SBATCH --output=../../output/07_ptwas/06_formatted_ld_panel/sbatch_output/array_%a.out
#SBATCH --error=../../output/07_ptwas/06_formatted_ld_panel/sbatch_output/array_%a.err
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --partition=nih_s10
#SBATCH --account=smontgom
#SBATCH --time=0-12:00:00

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate colocalization

python3 06_format_ld_panel.py $SLURM_ARRAY_TASK_ID
