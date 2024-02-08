#!/usr/bin/env bash
#SBATCH -J 402a_aggregation_sample
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

run_name="CCI_CellClass_L1_conf_malign"

input_file="${work_dir}/output/${run_name}/401_combine_samples/401_samples_interactions_mvoted.rds"
metadata="${work_dir}/output/${run_name}/000_data/bw_gbm_regional_study__metadata.rds"
min_cells=100
min_frac_samples=0.5
annot="CCI_CellClass_L1"
sample_varname="Sample"
patient_varname="Patient"
condition_varname="Region_Grouped"

output_dir="${work_dir}/output/${run_name}/402_aggregation"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"
Rscript "${work_dir}/scripts/402a_aggregation_sample.R" \
    --output_dir ${output_dir} \
    --input_file ${input_file} \
    --metadata ${metadata} \
    --min_cells ${min_cells} \
    --min_frac_samples ${min_frac_samples} \
    --annot ${annot} \
    --condition_varname ${condition_varname} \
    --sample_varname ${sample_varname}

echo "COMPLETED!"