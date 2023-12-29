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
# module load apptainer

# Local
base_dir="/cluster/projects/gaitigroup/Users/Joan/"
nf_exec="${HOME}/nextflow-23.04.3-all"
work_dir="${base_dir}/nf_work_cci"
nf_profile="slurm"

echo "Create work directory if not existing..."
mkdir -p $work_dir

project_dir="${base_dir}/scrnaseq-cellcomm"
outdir="nf-logs"

# TODO: if new run, then remove following
# rm -rf "${project_dir}/.nextflow"

echo "Setting user parameters..."
# input_file="/cluster/projects/gaitigroup/Users/Joan/002_Project_GBM/001_data/gbm_regional_study.rds"
# input_file=${project_dir}/output/CellClass_L4_min3_types_rerun/000_data/seurat_annot_adapted.rds
# input_file="${project_dir}/output/CellClass_L4_min3_types_rerun/000_data/seurat_annot_adapted_reduced_size.rds"
input_file="${project_dir}/001_data/test_data/example_data.rds"
# sample_dir="${project_dir}/output/CellClass_L4_min3_types_rerun/000_data/split_by_Sample"
# sample_dir="${project_dir}/output/CellClass_L4_min3_types_rerun/100_preprocessing"
samples_oi="/cluster/projects/gaitigroup/Users/Joan/002_Project_GBM/000_misc/samples_oi.txt"

metadata_csv="${project_dir}/output/CellClass_L4_min3_types_rerun/000_data/seurat_annot_adapted__metadata.csv"
metadata_rds="${project_dir}/output/CellClass_L4_min3_types_rerun/000_data/seurat_annot_adapted__metadata.rds"

skip_add_annot=1
skip_reduction=0
skip_preprocessing=0

split_varname="Sample"
annot="CellClass_L2"
min_cells=100
min_cell_types=3
n_perm=5
min_pct=0.10
alpha=0.05

output_run_name="test_pipeline"

# Databases of interactions
interactions_db="${project_dir}/001_data/interactions_db_v2"
cellphone_db="${interactions_db}/cellphonedb_custom/cellphonedb_12_18_2023_120229.zip"
cellchat_db="${interactions_db}/cellchat_db.rds"
liana_db="${interactions_db}/liana_db.rds"
liana_db_csv="${interactions_db}/liana_db.csv"
ref_db="${interactions_db}/ref_db.rds"
# interactions_db="${interactions_db}/001_data/interactions_db/interactions_ref.rds"

# Create output directory if not existing
mkdir -p "${project_dir}/output/${output_run_name}"
mkdir -p "${outdir}"

# CellClass_L4_min3_types: -resume "e9fdd3d9-27c0-4eae-a15e-1950ddb86eb2" \
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
    --output_run_name ${output_run_name} \
    --cellphone_db ${cellphone_db} \
    --cellchat_db ${cellchat_db} \
    --liana_db ${liana_db} \
    --liana_db_csv ${liana_db_csv} \
    --ref_db $ref_db \
    --outdir ${outdir} \
    --alpha $alpha #\
    # --skip_add_annot
    # --metadata_csv $metadata_csv \
    # --metadata_rds $metadata_rds \
echo "Done!"
