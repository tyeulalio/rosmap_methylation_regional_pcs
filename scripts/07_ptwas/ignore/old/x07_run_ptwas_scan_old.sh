#!/usr/bin/env bash

#SBATCH --job-name=ptwas_scan
#SBATCH --output=../../output/07_ptwas/07_ptwas_scan/array_%a.out
#SBATCH --error=../../output/07_ptwas/07_ptwas_scan/array_%a.err
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --partition=nih_s10
#SBATCH --account=smontgom
#SBATCH --time=7-00:00:00

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate colocalization

./07_ptwas_scan.sh
