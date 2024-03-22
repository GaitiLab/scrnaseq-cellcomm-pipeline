#!/usr/bin/env bash
#SBATCH -J launch_cci_pipeline_CCI_CellClass_L1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

module load java/18
# module load apptainer

# Local
base_dir="/cluster/projects/gaitigroup/Users/Joan/"
nf_exec="${HOME}/nextflow-23.04.3-all"
work_dir="${base_dir}/nf_work_cci_L1"
nf_profile="slurm"

echo "Create work directory if not existing..."
mkdir -p $work_dir

project_dir="${base_dir}/scrnaseq-cellcomm"

echo "PIPELINE CONFIGURATION..."
# General
run_name="CCI_CellClass_L1_w_agg"
approach=6

# Inputs 
input_file="/cluster/projects/gaitigroup/Data/GBM/processed_data/gbm_regional_study.rds"

# Pre-processing
split_varname="Sample"
annot="CCI_CellClass_L1"
condition_varname="Region_Grouped"
patient_varname="Patient"
min_patients=2
min_cells=100
min_cell_types=2

# Cell-cell interactions
n_perm=1000
min_pct=0.10
alpha=0.05

# Post-processing/formatting
meta_vars_oi="${project_dir}/000_misc/meta_vars_oi.txt"

# Databases of interactions
interactions_db="${project_dir}/data/interactions_db"
cellphone_db="${interactions_db}/cellphonedb.zip"
cellchat_db="${interactions_db}/cellchat_db.rds"
liana_db="${interactions_db}/liana_db.rds"
liana_db_csv="${interactions_db}/cell2cell_db.csv"
ref_db="${interactions_db}/ref_db.rds"

# Create output directory if not existing
mkdir -p "${project_dir}/output/${run_name}"
mkdir -p "${project_dir}/nf-logs"

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
    --run_name ${run_name} \
    --cellphone_db ${cellphone_db} \
    --cellchat_db ${cellchat_db} \
    --liana_db ${liana_db} \
    --liana_db_csv ${liana_db_csv} \
    --ref_db $ref_db \
    --alpha $alpha \
    --meta_vars_oi $meta_vars_oi \
    --approach $approach \
    --condition_varname $condition_varname \
    --aggregate_patients \
    --patient_varname $patient_varname \
    --min_patients $min_patients \
    --skip_downsampling
echo "Done!"
