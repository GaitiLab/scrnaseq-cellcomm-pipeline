#!/usr/bin/env bash

#SBATCH -J 201_cci_liana
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "standard_env"

n_perm=100
ident_col="CellClass_L4"
n_cells=100

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm

resource="${work_dir}/001_data/interactions_db_v2/liana_updated_test.rds"

sample="${work_dir}/output/CellClass_L4_min3_types/100_preprocessing/seurat/6419_cortex.rds"
output_dir="${work_dir}/output/test/"

Rscript "$work_dir/scripts/201_cci_liana.R" \
    --output_dir ${output_dir} \
    --ident_col $ident_col \
    --resource $resource \
    --n_perm $n_perm \
    --gene_expr $sample \
    --n_cells $n_cells