#!/usr/bin/env bash
#SBATCH -J extract_real_paths
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

run_name="CellClass_L4_min3_types"
project_dir="/cluster/projects/gaitigroup/Users/Joan/h4h-cell-cell-interactions/output/${run_name}/"

# # Data - preprocessing mtx
# data_dir="/cluster/projects/gaitigroup/Users/Joan/002_Project_GBM/001_data/preprocessing/"
# for folder in $data_dir/*; do
#     for subfolder in $folder/*; do
#         echo $(realpath ${subfolder})
#     done
# done >"${project_dir}/preprocessing.txt"

# # Cell-Cell-interactions
# for folder in $project_dir/*; do
#     for subfolder in $folder/*; do
#         echo $(realpath ${subfolder})
#     done
# done > "${project_dir}/output_files.txt"

global_project_dir="/cluster/projects/gaitigroup/Users/Joan/002_Project_GBM/001_data/preprocessing"
for folder in $global_project_dir/*; do
    for subfolder in $folder/*; do
        if [[ $subfolder != $(realpath ${subfolder}) ]]; then
            echo $(realpath ${subfolder})
        fi
    done
done > "${project_dir}/preprocessing.txt"