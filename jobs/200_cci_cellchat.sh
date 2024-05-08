#!/usr/bin/env bash
#SBATCH -J 200_cci_cellchat_6234_2895153_B
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=veryhimem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70G
#SBATCH --time=1-12:00:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

n_perm=1000
ident_col="CCI_CellClass_L2_2"

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm

resource="${work_dir}/data/interactions_db/cellchat_db.rds"
sample="${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/100_preprocessing/seurat/6509_cortex.rds"
output_dir="${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/200_cci_cellchat"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

Rscript "$work_dir/scripts/200_cci_cellchat.R" \
    --output_dir ${output_dir} \
    --ident_col $ident_col \
    --resource $resource \
    --n_perm $n_perm \
    --gene_expr $sample \
    --n_cores $SLURM_CPUS_PER_TASK
EOFs