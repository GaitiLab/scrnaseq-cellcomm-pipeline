#!/usr/bin/env bash
#SBATCH -J 501_heatmap
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=${base_dir}/Joan/scrnaseq-cellcomm

condition_varname="Region_Grouped"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

for annot in "CCI_CellClass_L1" "CCI_CellClass_L2"
do 
    run_name="${annot}_conf_min50"

    # Set output folder
    output_dir="${work_dir}/output/${run_name}/501_heatmap_n_interactions"
    for agg_level in "sample" "patient"
    do 
        input_file="${work_dir}/output/${run_name}/402_aggregation/402_${agg_level}_interactions_mvoted_w_filters.rds"

        Rscript "$work_dir/scripts/501_heatmap.R" \
            --input_file $input_file \
            --annot $annot \
            --output_dir $output_dir \
            --condition_varname $condition_varname \
            --agg_level $agg_level 
    done
done 
echo "COMPLETED!"