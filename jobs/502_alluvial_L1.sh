#!/usr/bin/env bash
#SBATCH -J 502_alluvial_L1
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

interactions_db="${work_dir}/data/interactions_db/ref_db.rds"
condition_varname="Region_Grouped"

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

annot="CCI_CellClass_L1"

# Set output folder
run_name="${annot}_conf_min50"
output_dir="${work_dir}/output/${run_name}/502_alluvial"
for agg_level in "sample" "patient"
do 
    input_file="${work_dir}/output/${run_name}/402_aggregation/402_${agg_level}_interactions_mvoted_w_filters.rds"

    Rscript "$work_dir/scripts/502_alluvial_L1.R" \
        --input_file $input_file \
        --annot $annot \
        --output_dir $output_dir \
        --condition_varname $condition_varname \
        --agg_level $agg_level \
        --interactions_db $interactions_db
done

echo "COMPLETED!"