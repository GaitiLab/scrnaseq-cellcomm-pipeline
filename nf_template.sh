#!/usr/bin/env bash
#SBATCH -J launch_cci_pipeline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
##SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

module load java/18


base_dir="/cluster/projects/gaitigroup/Users/Joan/"
# base_dir="/Users/joankant/Desktop/gaitigroup/Users/Joan"
project_dir="${base_dir}/scrnaseq-cellcomm"

# ---- PIPELINE CONFIGURATION ---- #
# Input seurat file
input_file="${project_dir}/data/example_data.rds"

# Output directory
output_dir="${project_dir}/test_pipeline"

init_step=1

# Pre-processing
sample_var="Sample"
annot="seurat_annotations"
condition_var="Condition"
patient_var="Patient"
min_patients=2
min_cells=70
is_confident=0

# Cell-cell interactions
n_perm=10
min_pct=0.10
alpha=0.05


# ----  NEXTFLOW CONFIGURATION ---- #
# Path to nextflow executable
nf_exec="${HOME}/nextflow-23.04.3-all"

# Work directory - all executed tasks (processes) are stored here
work_dir="${project_dir}/nf-work"
nf_profile="conda"
# Output directory for: trace, report + timeline by NextFlow
outdir="${project_dir}/nf-logs"

# Create directories
mkdir -p "${output_dir}"
mkdir -p "${outdir}"

echo "Running pipeline..."
# # Start the pipeline
# ${nf_exec} run ${project_dir} -with-report -with-trace \
${nf_exec} run ${project_dir} \
    -profile ${nf_profile} \
    -w ${work_dir} \
    --input_file $input_file \
    --sample_var ${sample_var} \
    --annot ${annot} \
    --min_cells ${min_cells} \
    --n_perm ${n_perm} \
    --min_pct ${min_pct} \
    --alpha $alpha \
    --init_step $init_step \
    --condition_var $condition_var \
    --patient_var $patient_var \
    --min_patients $min_patients \
    --is_confident ${is_confident} \
    --outdir ${outdir} \
    --output_dir ${output_dir}

echo "Done!"
