#!/usr/bin/env bash

#SBATCH --job-name=ptwas_estimate
#SBATCH --output=/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/09_ptwas_est/sbatch_output/array_%a.out
#SBATCH --error=/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/09_ptwas_est/sbatch_output/array_%a.err
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --partition=nih_s10
#SBATCH --account=smontgom
#SBATCH --time=0-12:00:00
#SBATCH --wait

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate r 


# arguments:
# 1: summary_types = ['avgs', 'pcs', 'cpgs']
# 2: region_types = ['full_gene', 'genome', 'astro_enh', 'promoters']
# 3: cell_types = ['astro', 'endo', 'neuron', 'oligo_opc', 'bulk']

# don't change these for now
summary_type="pcs"
region_type="full_gene"

#SLURM_ARRAY_TASK_ID=0

# select summary type here using the job array number
# remember this is zero-indexed so 0=astro, ..., 4=bulk
cell_types=("bulk" "astro" "endo" "neuron" "oligo_opc")
cell_type=${cell_types[$SLURM_ARRAY_TASK_ID]}

# process all cell_types for one region type, summary type
Rscript 09_ptwas_est.R $summary_type $region_type $cell_type
