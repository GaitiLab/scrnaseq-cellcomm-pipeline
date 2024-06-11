#!/usr/bin/env bash
#SBATCH -J launch_cci_pipeline_LP_IMM_10x
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiaoyi.chen@uhn.ca
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out
module load java/18

base_dir="/cluster/projects/gaitigroup/Users/Joan/"
project_dir="${base_dir}/scrnaseq-cellcomm-pipeline"

# ---- PIPELINE CONFIGURATION ---- #
# Input seurat file
input_file=/cluster/projects/gaitigroup/Users/Jiaoyi/analysis/Annotations/LP_IMM_3B13mut_annotation_scvi.rds

# Output directory
output_dir="${project_dir}/output/LP_IMM_perSample"

init_step=1

# Pre-processing
sample_var="Sample"
annot="CellClass_L3_LP"
condition_var="Mutation"
patient_var="Sample"
min_patients=2
min_cells=5

# Cell-cell interactions
n_perm=1000
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
${nf_exec} run ${project_dir} -resume \
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
    --outdir ${outdir} \
    --output_dir ${output_dir}

echo "Done!"
