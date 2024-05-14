#!/usr/bin/env bash
#SBATCH -J 401_combine_samples
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

output_dir="${work_dir}/output_Jiaoyi/401_combine_samples"
input_dir="${work_dir}/output_Jiaoyi/400_consensus"
metadata="/cluster/projects/gaitigroup/Users/Jiaoyi/scrnaseq-cellcomm/output/cci_scvi_merged_annotation_perSample_merged_CellClassL1_Apr12/000_data/merged_3B13mut_annotation_scvi_V2__metadata.rds"
# meta_vars_oi="${work_dir}/000_misc/meta_vars_oi.txt"
sample_varname="Sample"
patient_varname="Sample"
condition_varname="Mutation"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

Rscript "${work_dir}/scripts/401_combine_samples.R" \
    --output_dir $output_dir \
    --input_dir ${input_dir} \
    --metadata  ${metadata} \
    --condition_varname ${condition_varname} \
    --sample_varname ${sample_varname} \
    --patient_varname ${patient_varname}