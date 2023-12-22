#!/usr/bin/env bash
#SBATCH -J launch_200_cci_cellchat
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

job_min=1

n_perm=1000
ident_col="CellClass_L4"
n_cores=8

base_dir="${HOME}/Desktop/gaitigroup/Users"
# base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/h4h-cell-cell-interactions

resource="${work_dir}/001_data_local/interactions_db_v2/cellchat_db.rds"

sample_dir="${work_dir}/output/CellClass_L4_min3_types_rerun/100_preprocessing/seurat"
output_dir="${work_dir}/output/CellClass_L4_min3_types_rerun/200_cci_cellchat/"

sample_ids="${work_dir}/000_misc_local/rerun_cellchat.txt"
# # Determine job array limits
# # A. Determine number of files
# # job_max=$(ls -d -- $sample_dir/* | wc -l) 2>/dev/null
# # B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
# # job_max=1

# echo $job_max

# sbatch <<EOF
# #!/usr/bin/env bash

# #SBATCH -J 200_cci_cellchat
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=joan.kant@uhn.ca
# #SBATCH --partition=veryhimem
# #SBATCH --ntasks=1
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=${n_cores}
# #SBATCH --mem=64G
# #SBATCH --time=16:00:00
# #SBATCH --output=slurm_out/%x_%A_%a.out
# #SBATCH --error=slurm_out/%x_%A_%a.out
# #SBATCH --array=${job_min}-${job_max}

# echo "Activating conda environment..."
# source "\$HOME/miniforge3/bin/activate" "standard_env"

# SLURM_ARRAY_TASK_ID=1
# sample_id=\$(sed -n \${SLURM_ARRAY_TASK_ID}p ${sample_ids})
# sample="${sample_dir}/\${sample_id}"

# echo "Sample ID: \${sample_id}"
# echo "File:      \${sample}"
# Rscript "$work_dir/scripts/200_cci_cellchat.R" \
#     --output_dir ${output_dir} \
#     --ident_col $ident_col \
#     --resource $resource \
#     --n_perm $n_perm \
#     --gene_expr \${sample} \
#     --n_cores $n_cores
# EOF


SLURM_ARRAY_TASK_ID=1
sample_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${sample_ids})
sample="${sample_dir}/${sample_id}"

echo $sample > test.txt
sample2="/Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/output/CellClass_L4_min3_types_rerun/100_preprocessing/seurat/6234_2895153_B.rds"
echo $sample > test2.txt

if [ $sample == $sample2 ]; then 
    echo "same"
else
    echo "different"
fi


diff  -w <(echo "$sample"  ) <(echo "$sample2")

# Rscript "$work_dir/scripts/200_cci_cellchat.R" \
#     --output_dir ${output_dir} \
#     --ident_col $ident_col \
#     --resource $resource \
#     --n_perm $n_perm \
#     --gene_expr "${sample}" \
#     --n_cores $n_cores