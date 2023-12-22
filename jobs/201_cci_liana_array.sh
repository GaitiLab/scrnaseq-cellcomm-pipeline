#!/usr/bin/env bash
#SBATCH -J 201_cci_liana
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
ident_col="CellClass_L2"
n_cells=200

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/h4h-cell-cell-interactions

resource="${work_dir}/data/interactions_db/custom_liana.rds"

# /Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/final_output/CellClass_L2_all/100_preprocessing/seurat

# /Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/final_output/CellClass_L2_all/100_preprocessing/seurat
sample_dir="${work_dir}/final_output/${ident_col}_all/100_preprocessing/seurat"
output_dir="${work_dir}/final_output/${ident_col}_all/201_cci_liana/"

# Determine job array limits
# A. Determine number of files
job_max=$(ls -d -- $sample_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
# job_max=1

echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 201_cci_liana
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "standard_env"

sample=\$(ls -d -- $sample_dir/* | sed -n \${SLURM_ARRAY_TASK_ID}p)

Rscript "$work_dir/scripts/201_cci_liana.R" \
    --output_dir ${output_dir} \
    --ident_col $ident_col \
    --resource $resource \
    --n_perm $n_perm \
    --gene_expr \$sample \
    --n_cells $n_cells
EOF