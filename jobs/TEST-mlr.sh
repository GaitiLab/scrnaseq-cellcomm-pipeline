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
    output_dir="${work_dir}/output/${annot}/TEST-mlr"
    
    for type_of_voting in "lenient_voting" "stringent_voting"
    do 
        interactions_n_cells="${work_dir}/output/${annot}/TEST-ncells_impact/all_n_interactions_by_sample_${type_of_voting}.rds"
        potential_interactions="${work_dir}/output/${annot}/TEST-potential_vs_detected/all_n_interactions_by_sample_${type_of_voting}.rds"
        genes_per_class="${work_dir}/output/${annot}/TEST-potential_vs_detected/genes_per_class_${type_of_voting}.rds"

        Rscript "$work_dir/scripts/TEST-mlr.R" \
            --annot $annot \
            --output_dir $output_dir \
            --interactions_db $interactions_db \
            --meta $meta \
            --type_of_voting $type_of_voting \
            --interactions_n_cells $interactions_n_cells \
            --potential_interactions $potential_interactions \
            --genes_per_class $genes_per_class

            
    done
done 
echo "COMPLETED!"