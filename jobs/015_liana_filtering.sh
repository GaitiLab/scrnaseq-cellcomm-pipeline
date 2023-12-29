#!/usr/bin/env bash

#SBATCH -J 015_liana_filtering
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

source "${HOME}/miniforge3/bin/activate" "cci"

work_dir="/cluster/projects/gaitigroup/Users/Joan/scrnaseq-cellcomm"
output_dir="${work_dir}/output"
ref_dir="${work_dir}/000_misc/references"

Rscript "${work_dir}/scripts/015_liana_filtering.R" \
    --output_dir ${output_dir} \
    --ref_dir ${ref_dir}
