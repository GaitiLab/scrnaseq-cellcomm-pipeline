#!/usr/bin/env bash

#SBATCH -J compute_avg_expr_and_pseudobulk
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=01:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

source "${HOME}/miniforge3/bin/activate" "cci"

work_dir="/cluster/projects/gaitigroup/Users/Joan/scrnaseq-cellcomm"
input_file="/cluster/projects/gaitigroup/Data/GBM/processed_data/gbm_regional_study.rds"
output_dir="${work_dir}/output/average_expression_and_pseudobulk"

Rscript "${work_dir}/scripts/UTILS-compute_avg_expr_and_pseudobulk.R" \
    --input_file $input_file \
    --output_dir ${output_dir}