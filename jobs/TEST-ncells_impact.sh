#!/usr/bin/env bash
base_dir="${HOME}/Desktop/gaitigroup/Users"
# base_dir="/cluster/projects/gaitigroup/Users"
work_dir=${base_dir}/Joan/scrnaseq-cellcomm

interactions_db="${work_dir}/data/interactions_db/ref_db.rds"
meta="${work_dir}/output/CCI_CellClass_L2/000_data/bw_gbm_regional_study__metadata.rds"
percentile=0.9

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"
for annot in "CCI_CellClass_L1" "CCI_CellClass_L2"
do 
    # Set output folder
    output_dir="${work_dir}/output/${annot}/TEST-ncells_impact"
    interactions="${work_dir}/output/${annot}/401_combine_samples/401_samples_interactions_mvoted.rds"

        Rscript "$work_dir/scripts/TEST-ncells_impact.R" \
            --interactions $interactions \
            --annot $annot \
            --output_dir $output_dir \
            --meta $meta \
            --percentile $percentile
done 
echo "COMPLETED!"