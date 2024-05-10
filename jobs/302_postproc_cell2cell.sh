#!/usr/bin/env bash
#SBATCH -J launch_302_postproc_cell2cell
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


output_dir="${work_dir}/output_Jiaoyi/302_postproc_cell2cell"
ref_db="${work_dir}/001_data/interactions_db_v2/ref_db.rds"
sample_dir="/cluster/projects/gaitigroup/Users/Jiaoyi/scrnaseq-cellcomm/output/cci_scvi_merged_annotation_perSample_merged_CellClassL1_Apr12/202_cci_cell2cell"

# etermine job array limits
job_max=$(find $sample_dir -type f -name 'cell2cell__*.csv' | wc -l) 2>/dev/null

echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 302_postproc_cell2cell
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

input_file=\$(find $sample_dir -type f -name 'cell2cell__*.csv' | sed -n \${SLURM_ARRAY_TASK_ID}p)

base_name=\$(basename \$input_file)
tmp=\${base_name#"cell2cell__"}
sample_id=\${tmp%".csv"}

Rscript "$work_dir/scripts/302_postproc_cell2cell.R" \
    --output_dir ${output_dir} \
    --input_interactions \${input_file} \
    --ref_db ${ref_db} \
    --sample_id \${sample_id}
EOF