#!/usr/bin/env bash
#SBATCH -J launch_303_postproc_cpdb
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
work_dir=$base_dir/Joan/scrnaseq-cellcomm


output_dir="${work_dir}/output_Jiaoyi/303_postproc_cpdb"
ref_db="${work_dir}/data/interactions_db/ref_db.rds"
sample_dir="/cluster/projects/gaitigroup/Users/Jiaoyi/scrnaseq-cellcomm/output/cci_scvi_merged_annotation_perSample_merged_CellClassL1_Apr12/203_cci_cpdb"

# etermine job array limits
job_max=$(find $sample_dir -type f -name '*_metadata.tsv' | wc -l) 2>/dev/null

echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 303_postproc_cpdb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "cci"

input_file=\$(find $sample_dir -type f -name '*_metadata.tsv' | sed -n \${SLURM_ARRAY_TASK_ID}p)

base_name=\$(basename \$input_file)
sample_id=\${base_name%"_metadata.tsv"}

# Required files
interaction_scores="${sample_dir}/statistical_analysis_interaction_scores__\${sample_id}.txt"
pval="${sample_dir}/statistical_analysis_pvalues__\${sample_id}.txt"
sign_means="${sample_dir}/statistical_analysis_significant_means__\${sample_id}.txt"
means="${sample_dir}/statistical_analysis_means__\${sample_id}.txt"

Rscript "$work_dir/scripts/303_postproc_cellphonedb.R" \
    --output_dir ${output_dir} \
    --ref_db ${ref_db} \
    --sample_id \${sample_id} \
    --interaction_scores \${interaction_scores} \
    --pval \${pval} \
    --sign_means \${sign_means} \
    --means \${means}
EOF