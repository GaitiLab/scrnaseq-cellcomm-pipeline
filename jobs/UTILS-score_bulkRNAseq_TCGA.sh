#!/usr/bin/env bash
#SBATCH -J UTILS-score_bulkRNAseq_TCGA
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

source "${HOME}/miniforge3/bin/activate" "cci"

work_dir="/cluster/projects/gaitigroup/Users/Joan/scrnaseq-cellcomm"
output_dir="${work_dir}/output/score_bulkRNAseq"

meta="/cluster/projects/gaitigroup/Data/GBM/public_data/TCGA_bulkRNAseq/clinical.project-tcga-gbm.2024-04-17/clinical.tsv"
gene_exp="/cluster/projects/gaitigroup/Data/GBM/public_data/TCGA_bulkRNAseq/gdac.broadinstitute.org_GBM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/GBM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt"

suffix="CIBERSORT"
signatures="${work_dir}/000_misc/gene_lists/CIBERSORT_cell_type_signatures_20240429.rds"

Rscript "${work_dir}/scripts/UTILS-score_bulkRNAseq_TCGA.R" \
    --output_dir ${output_dir} \
    --signatures ${signatures} \
    --meta ${meta} \
    --gene_exp ${gene_exp} \
    --suffix ${suffix}


suffix="Neftel"
signatures="${work_dir}/000_misc/gene_lists/neftel_signatures.rds"
Rscript "${work_dir}/scripts/UTILS-score_bulkRNAseq_TCGA.R" \
    --output_dir ${output_dir} \
    --signatures ${signatures} \
    --meta ${meta} \
    --gene_exp ${gene_exp} \
    --suffix ${suffix}
