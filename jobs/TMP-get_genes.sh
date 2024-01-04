#!/usr/bin/env bash
#SBATCH -J tmp_get_genes
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm

input_file="${work_dir}/001_data/gbm_regional_study.rds"
output_dir="${work_dir}/output/"
source "$HOME/miniforge3/bin/activate" "standard_env"

Rscript "$work_dir/scripts/TMP-get_genes.R" \
    --input_file ${input_file} \
    --output_dir ${output_dir}

echo "COMPLETED!"