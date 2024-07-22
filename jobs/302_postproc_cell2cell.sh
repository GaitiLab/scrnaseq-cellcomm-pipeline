#!/usr/bin/env bash
#SBATCH -J launch_202_cci_cell2cell
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

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm-pipeline

run_dir="${base_dir}/Joan/GBM_CCI_Analysis/output/CCI_CellClass_L2_2_reassigned_samples_confident_only_FINAL"
sample_dir="${run_dir}/100_preprocessing/mtx"
output_dir="${run_dir}/302_postproc_cell2cell/"
input_dir="${run_dir}/202_cci_cell2cell"
ref_db="${work_dir}/data/interactions_db/ref_db.rds"

# Determine job array limits
# A. Determine number of files
job_max=$(ls -d -- $sample_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
# job_max=1
echo $job_max

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J 302_postproc_cell2cell
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Actsivating conda environment..."
source "\$HOME/miniforge3/bin/activate" "cci"

sample_id=\$(ls -d -- $sample_dir/* | sed -n \${SLURM_ARRAY_TASK_ID}p)
sample_id=\$(basename \${sample_id})

input_interactions_scores="${input_dir}/cell2cell__\${sample_id}__interaction_scores.csv"
input_interactions_pval="${input_dir}/cell2cell__\${sample_id}__pvalues.csv"

echo "Input scores: \${input_interactions_scores}"
echo "Input pvals: \${input_interactions_pval}"

echo "Sample ID: \${sample_id}"

Rscript "$work_dir/scripts/302_postproc_cell2cell.R" \
    --output_dir ${output_dir} \
    --input_interactions_scores "\${input_interactions_scores}" \
    --input_interactions_pval "\${input_interactions_pval}" \
    --ref_db ${ref_db} \
    --sample_id \${sample_id}
EOF