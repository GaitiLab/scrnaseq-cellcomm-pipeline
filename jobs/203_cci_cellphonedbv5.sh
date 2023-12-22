#!/usr/bin/env bash
#SBATCH -J launch_203_cci_cellphonedbv5
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

cpdb_file_path="${work_dir}/001_data/interactions_db_v2/cellphonedb_v5.0.0/cellphonedb.zip"
meta="${work_dir}/001_data/seurat_annot_adapted__metadata.csv"
annot="CellClass_L4"
nperm=1000
alpha=0.05
min_pct=0.1
nthreads=8

# TODO change to correct folders
sample_dir="${work_dir}/${annot}_all/100_preprocessing/mtx"
output_dir="${work_dir}/${annot}/203_cci_cellphonedbv5/"

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

#SBATCH -J 203_cci_cellphonedbv5
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
##SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${nthreads}
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Actsivating conda environment..."
source "\$HOME/miniforge3/bin/activate" "standard_env"

sample_dir=\$(ls -d -- $sample_dir/* | sed -n \${SLURM_ARRAY_TASK_ID}p)

python3 "$work_dir/Python/203_cci_cellphonedbv5.py" \
    --output_dir ${output_dir} \
    --meta $meta \
    --annot $annot \
    --cpdb_file_path $cpdb_file_path \
    --nperm $nperm \
    --input_dir \$sample_dir \
    --min_pct $min_pct \
    --threads \$SLURM_CPUS_PER_TASK \
EOF