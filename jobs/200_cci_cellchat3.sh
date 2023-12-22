#!/usr/bin/env bash
#SBATCH -J 200_cci_cellchat
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=veryhimem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out


n_perm=1000
ident_col="CellClass_L4"

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/h4h-cell-cell-interactions


resource="${work_dir}/001_data/interactions_db_v2/cellchat_db.rds"
# sample="${work_dir}/output/CellClass_L4_min3_types/100_preprocessing/seurat/6419_cortex.rds"
sample="${work_dir}/output/CellClass_L4_min3_types_rerun/100_preprocessing/seurat/6245_4972288_C.rds"
output_dir="${work_dir}/output/CellClass_L4_min3_types_rerun/200_cci_cellchat"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "standard_env"

Rscript "$work_dir/scripts/200_cci_cellchat.R" \
    --output_dir ${output_dir} \
    --ident_col $ident_col \
    --resource $resource \
    --n_perm $n_perm \
    --gene_expr $sample \
    --n_cores $SLURM_CPUS_PER_TASK
EOF