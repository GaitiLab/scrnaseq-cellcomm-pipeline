#!/usr/bin/env bash
#SBATCH -J launch_400a_consensus
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
run_name="CCI_CellClass_L2"

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm

run_dir="${work_dir}/output/${run_name}"
output_dir="${work_dir}/output/${run_name}/400_consensus"
alpha=0.05

sample_dir="${work_dir}/output/${run_name}/100_preprocessing/seurat"
# Determine job array limits
# A. Determine number of files
job_max=$(ls -d -- $sample_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
# job_max=1

echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 400a_consensus
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

Rscript "$work_dir/scripts/400a_consensus.R" \
    --run_dir ${run_dir} \
    --alpha ${alpha} \
    --sample_id \${sample_id} \
    --output_dir ${output_dir}
EOF