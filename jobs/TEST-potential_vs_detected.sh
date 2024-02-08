#!/usr/bin/env bash
base_dir="${HOME}/Desktop/gaitigroup/Users"
# base_dir="/cluster/projects/gaitigroup/Users"
work_dir=${base_dir}/Joan/scrnaseq-cellcomm

interactions_db="${work_dir}/data/interactions_db/ref_db.rds"
meta="${work_dir}/output/CCI_CellClass_L2/000_data/bw_gbm_regional_study__metadata.rds"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

for annot in "CCI_CellClass_L1" "CCI_CellClass_L2"
do 
    # Set output folder
    output_dir="${work_dir}/output/${annot}/TEST-potential_vs_detected"
    interactions="${work_dir}/output/${annot}/401_combine_samples/401_samples_interactions_mvoted.rds"
    gene_exp_dir="${work_dir}/output/${annot}/510_compute_avg_expr"

    for min_pct in 0 0.10
    do 

        Rscript "$work_dir/scripts/TEST-potential_vs_detected.R" \
            --interactions $interactions \
            --annot $annot \
            --output_dir $output_dir \
            --interactions_db $interactions_db \
            --meta $meta \
            --gene_exp_dir $gene_exp_dir \
            --min_pct $min_pct

            
    done
done 
echo "COMPLETED!"