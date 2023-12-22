#!/usr/bin/env bash
#SBATCH -J launch_cci_pipeline_CCI_CellClass_L2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=10-00:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

module load java/18
# module load apptainer

# Local
base_dir="/cluster/projects/gaitigroup/Users/Joan/"
nf_exec="${HOME}/nextflow-23.04.3-all"
work_dir="${base_dir}/nf_work_cci_L2"
nf_profile="slurm"

echo "Create work directory if not existing..."
mkdir -p $work_dir

project_dir="${base_dir}/h4h-cell-cell-interactions"

echo "PIPELINE CONFIGURATION..."
# General
output_run_name="CCI_CellClass_L2"
approach=3

# Inputs 
input_file="/cluster/projects/gaitigroup/Users/Joan/002_Project_GBM/001_data/gbm_regional_study.rds"

# Pre-processing
split_varname="Sample"
annot="CCI_CellClass_L2"
min_cells=100
min_cell_types=3
n_perm=1000
min_pct=0.10
alpha=0.05

# Databases of interactions
interactions_db="${project_dir}/001_data/interactions_db_v2"
cellphone_db="${interactions_db}/cellphonedb_custom/cellphonedb_12_18_2023_120229.zip"
cellchat_db="${interactions_db}/cellchat_db.rds"
liana_db="${interactions_db}/liana_db.rds"
liana_db_csv="${interactions_db}/liana_db.csv"
ref_db="${interactions_db}/ref_db.rds"

# Create output directory if not existing
mkdir -p "${project_dir}/output/${output_run_name}"

echo "Running pipeline..."
# # Start the pipeline
${nf_exec} run ${project_dir} -with-report -with-trace \
    -profile ${nf_profile} \
    -w ${work_dir} \
    --input_file $input_file \
    --split_varname ${split_varname} \
    --annot ${annot} \
    --min_cells ${min_cells} \
    --min_cell_types ${min_cell_types} \
    --n_perm ${n_perm} \
    --min_pct ${min_pct} \
    --output_run_name ${output_run_name} \
    --cellphone_db ${cellphone_db} \
    --cellchat_db ${cellchat_db} \
    --liana_db ${liana_db} \
    --liana_db_csv ${liana_db_csv} \
    --ref_db $ref_db \
    --alpha $alpha \
    --skip_add_annot
echo "Done!"
