#!/usr/bin/env bash
#SBATCH -J launch_401_combine_samples
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

run_name="CellClass_L4_min3_types"


# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/h4h-cell-cell-interactions

run_dir="${work_dir}/output/${run_name}/400_consensus"
output_dir="${work_dir}/output/${run_name}/401_combine_samples"
alpha=0.05
# celltype_oi="Neuron"
metadata="${work_dir}/001_data/seurat_annot_adapted__metadata.rds"

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 401_combine_samples
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

#input_file=\$(ls -d -- $interactions_dir/* | sed -n \${SLURM_ARRAY_TASK_ID}p)
# sample_id=\$(awk "NR==\${SLURM_ARRAY_TASK_ID}" ${sample_ids})
# input_file=$interactions_dir/\${sample_id}.rds

# is_stringent=\${stringency_options[\${SLURM_ARRAY_TASK_ID}]}

Rscript "$work_dir/scripts/401_combine_samples.R" \
    --input_dir ${run_dir} \
    --is_stringent \$SLURM_ARRAY_TASK_ID \
    --output_dir ${output_dir} \
    --metadata ${metadata}
EOF