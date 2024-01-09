#!/usr/bin/env bash

#SBATCH -J 001_reduce_seurat_object_size
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

source "${HOME}/miniforge3/bin/activate" "cci"

work_dir="/cluster/projects/gaitigroup/Users/Joan/scrnaseq-cellcomm"
input_file="${work_dir}/001_data/gbm_regional_study.rds"
output_dir="${work_dir}/output"


Rscript "${work_dir}/scripts/001_reduce_seurat_object_size.R" \
    --input_file $input_file \
    --output_dir ${output_dir}