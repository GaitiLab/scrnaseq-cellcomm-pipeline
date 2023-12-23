#!/usr/bin/env bash
#SBATCH -J run_script
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
work_dir=$base_dir/Joan/h4h-cell-cell-interactions

output_dir="${work_dir}/final_output/macrophage_tumor/"

input_dir="${work_dir}/final_output/macrophage_tumor/002_prep_sample/"
is_stringent=1

interactions="${work_dir}/final_output/macrophage_tumor/402_post_filtering/ccis_post_filtering__stringency_${is_stringent}_groupby_1.rds"

sample_ids="${work_dir}/final_output/macrophage_tumor/samples_oi_min200.txt"
# Determine job array limits
# A. Determine number of files
# job_max=$(ls -d -- $input_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
job_max=$(wc -l < "${sample_ids}")
# job_max=1

echo $job_max


sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 502_histograms
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:10:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "standard_env"

#input_file=\$(ls -d -- $input_dir/* | sed -n \${SLURM_ARRAY_TASK_ID}p)

sample_id=\$(awk "NR==\${SLURM_ARRAY_TASK_ID}" ${sample_ids})
filename="\${sample_id%.*}"


input_file=$input_dir/\${filename}__prepped.rds
echo \$input_file

Rscript "$work_dir/scripts/502_histograms.R" \
    --gene_expr \${input_file} \
    --output_dir ${output_dir} \
    --interactions ${interactions} \
    --is_stringent ${is_stringent}
EOF