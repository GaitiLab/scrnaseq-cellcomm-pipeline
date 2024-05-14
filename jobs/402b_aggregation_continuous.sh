#!/usr/bin/env bash
#SBATCH -J 402b_aggregation_continuous
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

input_file="${work_dir}/output_Jiaoyi/401_combine_samples//401_samples_interactions_agg_rank.rds"
condition_varname="Mutation"

output_dir="${work_dir}/output_Jiaoyi/402_aggregation"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

Rscript "${work_dir}/scripts/402b_aggregation_continuous.R" \
    --output_dir ${output_dir} \
    --input_file ${input_file} \
    --condition_varname ${condition_varname}

echo "COMPLETED!"