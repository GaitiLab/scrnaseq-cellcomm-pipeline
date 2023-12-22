#!/usr/bin/env bash

#SBATCH -J 000b_get_metadata
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

source "${HOME}/miniforge3/bin/activate" "standard_env"

work_dir="/cluster/projects/gaitigroup/Users/Joan/h4h-cell-cell-interactions"
input_file="/cluster/projects/gaitigroup/Users/Joan/002_Project_GBM/001_data/gbm_regional_study.rds"
output_dir="/cluster/projects/gaitigroup/Users/Joan/002_Project_GBM/001_data"

Rscript "${work_dir}/scripts/000b_get_metadata.R" \
    --input_file $input_file \
    --output_dir ${output_dir}