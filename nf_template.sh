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

run_name="test_pipeline"

# Setting up Nextflow
base_dir="/cluster/projects/gaitigroup/Users/Joan/"
nf_exec="${HOME}/nextflow-23.04.3-all"
work_dir="${base_dir}/nf_work_cci"
nf_profile="slurm"
# Output directory for: trace, report + timeline by NextFlow
outdir="nf-logs"
# Project directory
project_dir="${base_dir}/scrnaseq-cellcomm"

echo "Setting user parameters..."
input_file="${project_dir}/001_data/test_data/example_data.rds"

approach=5
split_varname="Sample"
annot="CellClass_L2"
min_cells=100
min_cell_types=3
n_perm=5
min_pct=0.10
alpha=0.05

# Databases of interactions (GitHub)
interactions_db="${project_dir}/data/interactions_db"
cellphone_db="${interactions_db}/cellphonedb_12_18_2023_120229.zip"
cellchat_db="${interactions_db}/cellchat_db.rds"
liana_db="${interactions_db}/liana_db.rds"
liana_db_csv="${interactions_db}/cell2cell_db.csv"
ref_db="${interactions_db}/ref_db.rds"

echo "Create work directory if not existing..."
mkdir -p $work_dir
mkdir -p $project_dir/$outdir
mkdir -p "${project_dir}/output/${run_name}"

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
    --outdir ${outdir} \
    --alpha $alpha 
echo "Done!"
