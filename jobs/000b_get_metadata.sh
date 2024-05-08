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

source "${HOME}/miniforge3/bin/activate" "cci"

work_dir="/cluster/projects/gaitigroup/Users/Joan/scrnaseq-cellcomm"
input_file="${work_dir}/001_data/gbm_regional_study.rds"
# input_file="/cluster/projects/gaitigroup/Data/GBM/public_data/gbmap_core.rds"
output_dir="${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/000_data"

Rscript "${work_dir}/scripts/000_get_metadata.R" \
    --input_file $input_file \
    --output_dir ${output_dir}