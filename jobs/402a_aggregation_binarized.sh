#!/usr/bin/env bash
#SBATCH -J 402a_aggregation_binarized
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

input_file="${work_dir}/output/${run_name}/401_combine_samples/401_samples_interactions_mvoted.rds"
metadata="${work_dir}/output/${run_name}/000_data/regional_study__metadata.rds"
annot="CCI_CellClass_L2_2"
condition_varname="Region"
min_patients=2

output_dir="${work_dir}/output/${run_name}/402_aggregation"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"
Rscript "${work_dir}/scripts/402a_aggregation_binarized.R" \
    --output_dir ${output_dir} \
    --input_file ${input_file} \
    --annot ${annot} \
    --condition_varname ${condition_varname} \
    --min_patients ${min_patients}

echo "COMPLETED!"