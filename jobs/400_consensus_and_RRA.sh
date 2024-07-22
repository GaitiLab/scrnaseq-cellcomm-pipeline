#!/usr/bin/env bash
#SBATCH -J launch_400_consensus_and_RRA.R
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

# run_dir="${work_dir}/output/multiome"
run_dir="${base_dir}/Joan/GBM_CCI_Analysis/output/CCI_CellClass_L2_2_reassigned_samples_confident_only_FINAL"

output_dir="${run_dir}/400_consensus_and_RRA"
alpha=0.05
n_perm=1000
sample_dir="${run_dir}/100_preprocessing/seurat"
# Determine job array limits
# A. Determine number of files
job_max=$(ls -d -- $sample_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
# job_max=1

echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 400_consensus_and_RRA.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "cci"

input_file=\$(ls -d -- $sample_dir/* | sed -n \${SLURM_ARRAY_TASK_ID}p)
filename=\$(basename -- "\$input_file")
sample_id="\${filename%.*}"

echo \$sample_id

Rscript "$work_dir/scripts/400_consensus_and_RRA.R" \
    --run_dir ${run_dir} \
    --alpha ${alpha} \
    --sample_id \${sample_id} \
    --output_dir ${output_dir} \
    --n_perm ${n_perm}
EOF