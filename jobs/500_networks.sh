#!/usr/bin/env bash
#SBATCH -J 500_networks
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
min_q=0.5
min_cells=100
has_loops=1
min_celltypes=2

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "cci"

for annot in "CCI_CellClass_L1" "CCI_CellClass_L2"
do 
    # Set output folder
    run_name="${annot}_conf_malign"
    output_dir="${work_dir}/output/${run_name}/500_networks"
    colors="${work_dir}/000_misc_local/${run_name}_network_colors.rds"
    meta="${work_dir}/output/${run_name}/000_data/bw_gbm_regional_study__metadata.rds"

    for agg_level in "sample" "patient"
    do 
        input_file="${work_dir}/output/${run_name}/402_aggregation/402_${agg_level}_interactions_mvoted_w_filters.rds"

        Rscript "$work_dir/scripts/500_networks.R" \
            --input_file $input_file \
            --annot $annot \
            --output_dir $output_dir \
            --agg_level $agg_level \
            --meta $meta \
            --interactions_db $interactions_db \
            --min_q $min_q \
            --min_cells $min_cells \
            --has_loops $has_loops \
            --min_celltypes $min_celltypes
    done
done 
echo "COMPLETED!"