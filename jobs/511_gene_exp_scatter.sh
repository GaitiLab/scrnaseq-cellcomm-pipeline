#!/usr/bin/env bash
#SBATCH -J 511_gene_exp_scatter
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=${base_dir}/Joan/scrnaseq-cellcomm

metadata="${work_dir}/output/CCI_CellClass_L2/000_data/bw_gbm_regional_study__metadata.rds"

interactions_db="${work_dir}/data/interactions_db/ref_db.rds"
min_pct_exp=10
gene_exp_dir="${work_dir}/output/${run_name}/510_compute_avg_expr"
interactions="${work_dir}/output/${run_name}/400_consensus/400_samples_interactions_mvoted_w_filters.rds"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

for annot in "CCI_CellClass_L1" "CCI_CellClass_L2"
do 
    run_name="${annot}_conf_malign"
    gene_exp_dir="${work_dir}/output/${run_name}/510_compute_avg_expr"
    for agg_level in "sample" "patient"
        do 
        echo "Annot=${annot}; agg_level=${agg_level}"
        interactions="${work_dir}/output/${run_name}/402_aggregation/402_${agg_level}_interactions_mvoted_w_filters.rds"
        output_dir="${work_dir}/output/${run_name}/511_gene_exp_scatter/${agg_level}"

        Rscript "$work_dir/scripts/511_gene_exp_scatter.R" \
            --meta "$metadata" \
            --output_dir "$output_dir" \
            --interactions_db "$interactions_db" \
            --min_pct_exp "$min_pct_exp" \
            --gene_exp_dir "$gene_exp_dir" \
            --interactions "$interactions"   \
            --annot $annot
    done
done
echo "COMPLETED!"