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

nperm=1000
annot="CellClass_L2"

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/h4h-cell-cell-interactions

resource="${work_dir}/data/interactions_db/custom_liana.csv"

meta="/cluster/projects/gaitigroup/Users/Joan/001_data/GBM/combined_metadata.csv"
# /Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/final_output/CellClass_L2_all/100_preprocessing/seurat

# /Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/final_output/CellClass_L2_all/100_preprocessing/seurat
sample_dir="${work_dir}/final_output/${annot}_all/100_preprocessing/mtx"
output_dir="${work_dir}/final_output/${annot}_all/202_cci_cell2cell/"

# Determine job array limits
# A. Determine number of files
job_max=$(ls -d -- $sample_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
# job_max=1
job_max=19
job_max=20
echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 202_cci_cell2cell
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=veryhimem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Actsivating conda environment..."
source "\$HOME/miniforge3/bin/activate" "standard_env"

sample_dir=\$(ls -d -- $sample_dir/* | sed -n \${SLURM_ARRAY_TASK_ID}p)

python3 "$work_dir/Python/202_cci_cell2cell.py" \
    --output_dir ${output_dir} \
    --meta $meta \
    --annot $annot \
    --interactions $resource \
    --nperm $nperm \
    --input_dir \$sample_dir
EOF