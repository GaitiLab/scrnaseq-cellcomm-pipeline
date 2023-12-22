#!/usr/bin/env bash

#SBATCH -J 000b_get_metadata
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

source "${HOME}/miniforge3/bin/activate" "standard_env"

work_dir="/cluster/projects/gaitigroup/Users/Joan/h4h-cell-cell-interactions"

Rscript ${work_dir}/TMP_scripts.R