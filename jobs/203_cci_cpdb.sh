#!/usr/bin/env bash
#SBATCH -J launch_203_cci_cpdb
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

base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm

cpdb_file_path="${work_dir}/data/interactions_db/cellphonedb.zip"
meta="${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/000_data/gbm_regional_study__metadata.csv"
annot="CCI_CellClass_L2_2"
nperm=1000
alpha=1.1
min_pct=0.1
nthreads=8

# TODO change to correct folders
# sample_dir="${work_dir}/${annot}_all/100_preprocessing/mtx"
output_dir="${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/203_cci_cpdb/"

# Determine job array limits
# A. Determine number of files
# job_max=$(ls -d -- $sample_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
job_max=1
# job_max=19
# job_max=20
sample_dir=${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/100_preprocessing/mtx/6509_cortex
echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 203_cci_cellphonedb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
##SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${nthreads}
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Actsivating conda environment..."
source "\$HOME/miniforge3/bin/activate" "cpdb"

sample_dir=\$(ls -d -- $sample_dir/* | sed -n \${SLURM_ARRAY_TASK_ID}p)

python3 "$work_dir/Python/203_cci_cpdb.py" \
    --output_dir ${output_dir} \
    --meta $meta \
    --annot $annot \
    --cpdb_file_path $cpdb_file_path \
    --n_perm $nperm \
    --input_dir $sample_dir \
    --min_pct $min_pct \
    --alpha ${alpha} \
    --threads \$SLURM_CPUS_PER_TASK
EOF