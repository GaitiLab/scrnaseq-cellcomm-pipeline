#!/usr/bin/env bash
#SBATCH -J launch_score_scRNAseq
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
work_dir=$base_dir/Joan/scrnaseq-cellcomm

signatures="${work_dir}/000_misc/gene_lists/neftel_signatures.rds"


# Our cohort
output_dir="${work_dir}/output/score_scRNAseq/our_cohort"
seurat_obj="${base_dir}/Data/GBM/processed_data/gbm_regional_study.rds"
samples_oi="${work_dir}/000_misc/gbm_regional_study__samples.txt"
sample_varname="Sample"
keep_confident=1
malignant_label="Malignant"
annot="CellClass_L1"

# Determine number of jobs
job_max=$(wc -l < "${samples_oi}")
echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J score_scRNAseq_our_cohort
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=01:00:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "cci"

sample_id=\$(awk "NR==\${SLURM_ARRAY_TASK_ID}" ${samples_oi})

echo Sample ID: \$sample_id
echo Input file path: \$input_file

Rscript ${work_dir}/scripts/UTILS-score_scRNAseq.R \
    --input_file \${input_file} \
    --output_dir ${output_dir} \
    --sample_id \${sample_id} \
    --sample_varname ${sample_varname} \
    --keep_confident ${keep_confident} \
    --malignant_label ${malignant_label} \
    --signatures ${signatures} \
    --annot ${annot}

EOF


# GBMap
output_dir="${work_dir}/output/score_scRNAseq/GBMap"
seurat_obj="${base_dir}/Data/GBM/public_data/gbmap_core.rds"
samples_oi="${work_dir}/000_misc/gbmap_core__samples.txt"
sample_varname="sample"
keep_confident=0
malignant_label="malignant cell"
annot="cell_type"

# Determine number of jobs
job_max=$(wc -l < "${samples_oi}")
echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J score_scRNAseq_GBMap
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=01:00:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "cci"

sample_id=\$(awk "NR==\${SLURM_ARRAY_TASK_ID}" ${samples_oi})

echo Sample ID: \$sample_id
echo Input file path: \$input_file

Rscript ${work_dir}/scripts/UTILS-score_scRNAseq.R \
    --input_file \${input_file} \
    --output_dir ${output_dir} \
    --sample_id \${sample_id} \
    --sample_varname ${sample_varname} \
    --keep_confident ${keep_confident} \
    --malignant_label ${malignant_label} \
    --signatures ${signatures} \
    --annot ${annot}

EOF