#!/usr/bin/env bash
#SBATCH -J GBM-create_figures
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

source "${HOME}/miniforge3/bin/activate" "cci"

work_dir="/cluster/projects/gaitigroup/Users/Joan/scrnaseq-cellcomm"
output_dir="/cluster/projects/gaitigroup/Users/Joan/002_Project_GBM/002_WIP_FIGURES/scRNAseq/CCI"


run_name="CCI_CellClass_L2_2_reassigned_samples_confident_only"
interactions_run_dir="${work_dir}/output/${run_name}"
condition_varname="Region"
unique_interactions="${work_dir}/output/${run_name}/${run_name}_unique_interactions_neuron_invasive_high.xlsx"
top_n=10
# source_oi="Neuron"
# target_oi="Invasive-high OPC_NPC1"
condition_oi="PT"
interactions_db="${work_dir}/data/interactions_db/ref_db.rds"
meta="${work_dir}/output/${run_name}/000_data/gbm_regional_study__metadata.rds"
avg_expr="${work_dir}/output/average_expression_and_pseudobulk/average_expression_by_Sample_CCI_CellClass_L2_2.rds"

Rscript "${work_dir}/scripts/GBM-create_figures.R" \
    --output_dir ${output_dir} \
    --interactions_run_dir $interactions_run_dir \
    --condition_varname $condition_varname \
    --unique_interactions $unique_interactions \
    --top_n $top_n \
    --condition_oi $condition_oi \
    --interactions_db $interactions_db \
    --meta $meta \
    --avg_expr $avg_expr