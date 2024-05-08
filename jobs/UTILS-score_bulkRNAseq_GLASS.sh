#!/usr/bin/env bash
#SBATCH -J UTILS-score_bulkRNAseq_GLASS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

source "${HOME}/miniforge3/bin/activate" "cci"

work_dir="/cluster/projects/gaitigroup/Users/Joan/scrnaseq-cellcomm"
output_dir="${work_dir}/output/score_bulkRNAseq"
meta="/cluster/projects/gaitigroup/Data/GBM/public_data/GLASS/analysis_rnaseq_pairs.csv"
gene_exp="/cluster/projects/gaitigroup/Data/GBM/public_data/GLASS/gene_tpm_matrix_all_samples.tsv"
status_info="/cluster/projects/gaitigroup/Data/GBM/public_data/GLASS/clinical_surgerie.csv"

suffix="CIBERSORT"
signatures="${work_dir}/000_misc/gene_lists/CIBERSORT_cell_type_signatures_20240429.rds"

Rscript "${work_dir}/scripts/UTILS-score_bulkRNAseq_GLASS.R" \
    --output_dir ${output_dir} \
    --signatures ${signatures} \
    --meta ${meta} \
    --gene_exp ${gene_exp} \
    --status_info ${status_info} \
    --suffix ${suffix}

suffix="Neftel"
signatures="${work_dir}/000_misc/gene_lists/neftel_signatures.rds"
Rscript "${work_dir}/scripts/UTILS-score_bulkRNAseq_GLASS.R" \
    --output_dir ${output_dir} \
    --signatures ${signatures} \
    --meta ${meta} \
    --gene_exp ${gene_exp} \
    --status_info ${status_info} \
    --suffix ${suffix}

