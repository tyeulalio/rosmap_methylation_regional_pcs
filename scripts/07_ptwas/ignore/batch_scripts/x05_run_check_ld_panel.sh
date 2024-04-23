#!/usr/bin/env bash

#SBATCH --job-name=check_ld_panel
#SBATCH --output=../../../output/07_ptwas/05_lifted_ld_panel/array_%a.out
#SBATCH --error=../../../output/07_ptwas/05_lifted_ld_panel/array_%a.err
#SBATCH --array=1-22
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
mamba activate base
mamba activate r
#source /home/eulalio/.bashrc
#source activate /home/eulalio/micromamba/envs/r


Rscript ../05_check_ld_panel.R $SLURM_ARRAY_TASK_ID
