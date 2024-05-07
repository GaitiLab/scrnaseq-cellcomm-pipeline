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

# Setting up Nextflow
base_dir="/cluster/projects/gaitigroup/Users/Joan/"
project_dir="${base_dir}/scrnaseq-cellcomm"

nf_exec="${HOME}/nextflow-23.04.3-all"
work_dir="${project_dir}/nf-work"
nf_profile="conda"
# Output directory for: trace, report + timeline by NextFlow
outdir="${project_dir}}/nf-logs"
# Project directory

echo "PIPELINE CONFIGURATION..."
# General
run_name="test_pipeline"
approach=6

# Inputs 
input_file="${project_dir}/data/example_data.rds"

# Downsampling
# downsampling_sheet="${project_dir}/output/downsampling_info.rds"
# num_repeats=3
# num_cells=2000

# Pre-processing
split_varname="Sample"
annot="seurat_annotations"
condition_varname="Condition"
patient_varname="Sample"
min_patients=2
min_cells=50
min_cell_types=2

# Cell-cell interactions
n_perm=100
min_pct=0.10
alpha=0.05

# Post-processing/formatting
# meta_vars_oi="${project_dir}/000_misc/meta_vars_oi.txt"

# Databases of interactions (GitHub)
interactions_db="${project_dir}/data/interactions_db"
cellphone_db="${interactions_db}/cellphonedb.zip"
cellphone_db="${interactions_db}/cellphonedb.zip"
cellchat_db="${interactions_db}/cellchat_db.rds"
liana_db="${interactions_db}/liana_db.rds"
liana_db_csv="${interactions_db}/cell2cell_db.csv"
ref_db="${interactions_db}/ref_db.rds"

echo "Create work directory if not existing..."
# Create output directory if not existing
mkdir -p "${project_dir}/output/${run_name}"
mkdir -p "${project_dir}/nf-logs"

echo "Running pipeline..."
# # Start the pipeline
${nf_exec} run ${project_dir} -with-report -with-trace -resume \
    -profile ${nf_profile} \
    -w ${work_dir} \
    --input_file $input_file \
    --split_varname ${split_varname} \
    --annot ${annot} \
    --min_cells ${min_cells} \
    --min_cell_types ${min_cell_types} \
    --n_perm ${n_perm} \
    --min_pct ${min_pct} \
    --run_name ${run_name} \
    --run_name ${run_name} \
    --cellphone_db ${cellphone_db} \
    --cellchat_db ${cellchat_db} \
    --liana_db ${liana_db} \
    --liana_db_csv ${liana_db_csv} \
    --ref_db $ref_db \
    --alpha $alpha \
    --approach $approach \
    --condition_varname $condition_varname \
    --patient_varname $patient_varname \
    --min_patients $min_patients \
    --num_cells ${num_cells} \
    --num_repeats ${num_repeats} \
    --skip_downsampling
echo "Done!"
