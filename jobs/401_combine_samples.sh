#!/usr/bin/env bash
#SBATCH -J 401_combine_samples
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir="$base_dir/Joan/scrnaseq-cellcomm-pipeline"

run_dir="${base_dir}/Joan/GBM_CCI_Analysis/output/CCI_CellClass_L2_2_reassigned_samples_confident_only_FINAL"

output_dir="${run_dir}/401_combine_samples"
input_dir="${run_dir}/400_consensus_and_RRA"
metadata="${run_dir}/000_data/gbm_regional_study__metadata.rds"
# meta_vars_oi="${work_dir}/000_misc/meta_vars_oi.txt"
sample_var="Sample"
patient_var="Patient"
condition_var="Region"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

Rscript "${work_dir}/scripts/401_combine_samples.R" \
    --output_dir $output_dir \
    --input_dir ${input_dir} \
    --metadata  ${metadata} \
    --condition_var ${condition_var} \
    --sample_var ${sample_var} \
    --patient_var ${patient_var}