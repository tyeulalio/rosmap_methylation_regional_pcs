#!/usr/bin/env bash

#SBATCH --job-name=sim_dmr
#SBATCH --output=../../output/03_simulation_dmrs/sbatch_output/array_%a.out
#SBATCH --error=../../output/03_simulation_dmrs/sbatch_output/array_%a.err
#SBATCH --array=1-160
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --partition=nih_s10
#SBATCH --account=smontgom
#SBATCH --time=0-3:00:00
#SBATCH --wait

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
#mamba activate R4.2
mamba activate r

conda info --envs

echo "Processing job " $SLURM_ARRAY_TASK_ID 

Rscript ./03_simulation_dmrs.R $SLURM_ARRAY_TASK_ID
