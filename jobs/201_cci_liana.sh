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
source "$HOME/miniforge3/bin/activate" "cci"

n_perm=1000
ident_col="CCI_CellClass_L2_2"
n_cells=50

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm

resource="${work_dir}/data/interactions_db/liana_db.rds"

sample="${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/100_preprocessing/seurat/6509_cortex.rds"
output_dir="${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/201_cci_liana/"

Rscript "$work_dir/scripts/201_cci_liana.R" \
    --output_dir ${output_dir} \
    --ident_col $ident_col \
    --resource $resource \
    --n_perm $n_perm \
    --gene_expr $sample \
    --n_cells $n_cells