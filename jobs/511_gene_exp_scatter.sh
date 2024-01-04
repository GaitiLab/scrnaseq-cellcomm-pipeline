#!/usr/bin/env bash
#SBATCH -J 511_gene_exp_scatter
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=${base_dir}/Joan/scrnaseq-cellcomm

run_name="CCI_CellClass_L2"

metadata="${work_dir}/output/${run_name}/000_data/gbm_regional_study__metadata.rds"

output_dir="${work_dir}/output/${run_name}/511_gene_exp_scatter"
interactions_db="${work_dir}/001_data/interactions_db_v2/ref_db.rds"
cutoff_quantile=0.90
min_pct_exp=10
res_threshold=1.96
gene_exp_dir="${work_dir}/output/${run_name}/510_compute_avg_expr"
interactions="${work_dir}/output/${run_name}/400_consensus/400_samples_interactions_mvoted_w_filters.rds"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

Rscript "$work_dir/scripts/511_gene_exp_scatter.R" \
    --meta "$metadata" \
    --output_dir "$output_dir" \
    --interactions_db "$interactions_db" \
    --cutoff_quantile "$cutoff_quantile" \
    --min_pct_exp "$min_pct_exp" \
    --res_threshold "$res_threshold" \
    --gene_exp_dir "$gene_exp_dir" \
    --interactions "$interactions"   

echo "COMPLETED!"