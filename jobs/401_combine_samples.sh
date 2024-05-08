#!/usr/bin/env bash
#SBATCH -J 400_combine_samples
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out


run_name="CCI_CellClass_L2_2_reassigned_samples_confident_only"

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm

output_dir="${work_dir}/output/${run_name}/401_combine_samples"
metadata="${work_dir}/output/${run_name}/000_data/gbm_regional_study__metadata.rds"
meta_vars_oi="${work_dir}/000_misc/meta_vars_oi.txt"
sample_varname="Sample"
patient_varname="Patient"
condition_varname="Region_Grouped"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

Rscript "${work_dir}/scripts/401_combine_samples.R" \
    --output_dir $output_dir \
    --input_dir ${work_dir}/output/${run_name}/400_consensus \
    --metadata  ${metadata} \
    --meta_vars_oi ${meta_vars_oi} \
    --condition_varname ${condition_varname} \
    --sample_varname ${sample_varname} \
    --patient_varname ${patient_varname}