#!/usr/bin/env bash
#SBATCH -J 402c_aggregation_integration
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=${base_dir}/Joan/scrnaseq-cellcomm

run_name="CCI_CellClass_L2_2_reassigned_samples_confident_only"

input_dir="${work_dir}/output_Jiaoyi/402_aggregation"
condition_varname="Mutation"

output_dir="${work_dir}/output_Jiaoyi/402_aggregation"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

Rscript "${work_dir}/scripts/402c_aggregation_integration.R" \
    --output_dir ${output_dir} \
    --input_dir ${input_dir} \
    --condition_varname ${condition_varname}

echo "COMPLETED!"