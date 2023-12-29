#!/usr/bin/env bash

#SBATCH -J 014_update_liana
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
cellchat_db="${work_dir}/001_data/interactions_db_v2/cellchat_liana_format.rds"
cpdb_db="${work_dir}/001_data/interactions_db_v2/cpdbv5_liana_format.rds"
liana_db="${work_dir}/001_data/interactions_db_v2/liana_db_base.rds"

Rscript "${work_dir}/scripts/014_update_liana.R" \
    --output_dir ${output_dir} \
    --cellchat_db ${cellchat_db} \
    --cpdb_db ${cpdb_db} \
    --liana_db ${liana_db}