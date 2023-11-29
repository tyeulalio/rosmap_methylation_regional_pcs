#!/usr/bin/env bash

#SBATCH --job-name=sim
#SBATCH --output=./sbatch_output/array_%a.out
#SBATCH --error=./sbatch_output/array_%a.err
#SBATCH --array=6-270
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --partition=nih_s10
#SBATCH --account=smontgom
#SBATCH --time=4-00:00:00

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate R4.2

conda info --envs

echo "Processing job " $SLURM_ARRAY_TASK_ID 

Rscript ./00_simulation_test.R $SLURM_ARRAY_TASK_ID
