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

run_dir="${work_dir}/final_output/CellClass_L2_all/401_combine_samples"
output_dir="${work_dir}/final_output/CellClass_L2_all/"
alpha=0.05
# celltype_oi="Neuron"

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 402_post_filtering
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
##SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=0-1

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "standard_env"

input_file="${run_dir}/401_samples_combined__stringency_\$SLURM_ARRAY_TASK_ID.rds"

Rscript "$work_dir/scripts/402_post_filtering.R" \
    --input_file \${input_file} \
    --is_stringent \$SLURM_ARRAY_TASK_ID \
    --output_dir ${output_dir}
EOF