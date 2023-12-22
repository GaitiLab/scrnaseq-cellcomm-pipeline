#!/usr/bin/env bash
#SBATCH -J run_script
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
##SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

job_min=1

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/h4h-cell-cell-interactions

output_dir="${work_dir}/final_output/CellClass_L2_all/"
input_dir="${work_dir}/final_output/CellClass_L2_all/100_preprocessing/seurat"

interactions_ref="${work_dir}/data/interactions_db/interactions_ref.rds"
# celltypes_oi="${work_dir}/data/celltypes_oi_neuron.txt"
interactions_oi="${work_dir}/final_output/CellClass_L2_all/402_post_filtering/ccis_post_filtering__stringency_1_groupby_1.rds"
min_pct_expressed=0
# Determine job array limits
# A. Determine number of files
job_max=$(ls -d -- $input_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
# job_max=1

echo $job_max


sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 501a_compute_avg_expr
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "standard_env"

input_file=\$(ls -d -- $input_dir/* | sed -n \${SLURM_ARRAY_TASK_ID}p)

# input_file="${run_dir}/401_samples_combined__stringency_\$SLURM_ARRAY_TASK_ID.rds"

Rscript "$work_dir/scripts/501a_compute_avg_expr.R" \
    --gene_expr \${input_file} \
    --interactions_ref ${interactions_ref} \
    --interactions_oi ${interactions_oi} \
    --min_pct_expressed ${min_pct_expressed} \
    --output_dir ${output_dir}
EOF