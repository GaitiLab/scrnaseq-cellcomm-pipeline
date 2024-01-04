#!/usr/bin/env bash

#SBATCH -J 013_get_cpdb
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
cpddb_dir="${work_dir}/000_misc/references/cellphonedb_v5.0.0"

Rscript "${work_dir}/scripts/013_get_cpdb.R" \
    --output_dir ${output_dir} \
    --cpddb_dir ${cpddb_dir}