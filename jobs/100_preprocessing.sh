#!/usr/bin/env bash
#SBATCH -J launch_preprocessing
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

annot="CellClass_L2"
min_cells=50
celltypes_oi=""
# first_n = 0, don't check for presence of specific cell tyeps
first_n=0

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm

interactions_db="${work_dir}/data/interactions_db/interactions_ref.rds"

sample_dir="/cluster/projects/gaitigroup/Users/Joan/001_data/GBM/split_by_Sample"
output_dir="${work_dir}/final_output/${annot}_all/100_preprocessing/"

# Determine job array limits
# A. Determine number of files
job_max=$(ls -d -- $sample_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
# job_max=1

echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 100_preprocessing
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

Rscript "$work_dir/scripts/100_preprocessing.R" \
    --input_file \${sample} \
    --output_dir ${output_dir} \
    --annot ${annot} \
    --min_cells ${min_cells} \
    --first_n ${first_n} \
    --interactions_db ${interactions_db}
EOF