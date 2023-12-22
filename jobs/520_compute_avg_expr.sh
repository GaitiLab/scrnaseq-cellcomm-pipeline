#!/usr/bin/env bash
#SBATCH -J launch_520_compute_avg_expr
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
n_cells=200

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/h4h-cell-cell-interactions

resource="${work_dir}/001_data/interactions_db/interactions_ref.rds"

sample_dir="${work_dir}/output/CellClass_L4_min3_types/100_preprocessing/seurat"
output_dir="${work_dir}/output/CellClass_L4_min3_types/520_compute_avg_expr/"

# Determine job array limits
# A. Determine number of files
job_max=$(ls -d -- $sample_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
# job_max=1

echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 520_compute_avg_expr
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
##SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "standard_env"

sample=\$(ls -d -- $sample_dir/* | sed -n \${SLURM_ARRAY_TASK_ID}p)

Rscript "$work_dir/scripts/520_compute_avg_expr.R" \
    --output_dir ${output_dir} \
    --interactions_db $resource \
    --gene_exp \$sample 
EOF