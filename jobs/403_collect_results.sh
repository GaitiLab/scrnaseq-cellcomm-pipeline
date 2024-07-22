#!/usr/bin/env bash
#SBATCH -J 403_collect_results
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

base_dir="${HOME}/Desktop/gaitigroup/Users"
# base_dir="/cluster/projects/gaitigroup/Users"
work_dir="$base_dir/Joan/scrnaseq-cellcomm-pipeline"

run_dir="${base_dir}/Joan/GBM_CCI_Analysis/output/CCI_CellClass_L2_2_reassigned_samples_confident_only_FINAL"

output_dir="${run_dir}/"
interactions_agg_integration="${run_dir}/402_aggregation_and_filtering/402c_filtering_aggregated_res.rds"
condition_var="Region"
output_name="interactions_summary"
alpha=0.05
echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

Rscript "${work_dir}/scripts/403_collect_results.R" \
    --output_dir $output_dir \
    --interactions_agg_integration ${interactions_agg_integration} \
    --condition_var ${condition_var} \
    --alpha ${alpha} \
    --output_name ${output_name}
