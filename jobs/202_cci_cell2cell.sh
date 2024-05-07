#!/usr/bin/env bash
#SBATCH -J 202_cci_cell2cell
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out


nperm=1000
annot="CCI_CellClass_L2_2"

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm
sample_id="6509_cortex"
resource="${work_dir}/data/interactions_db/cell2cell_db.csv"

meta="${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/000_data/gbm_regional_study__metadata.csv"
sample_dir="${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/100_preprocessing/mtx/6509_cortex"
output_dir="${work_dir}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/202_cci_cell2cell"

echo "Actsivating conda environment..."
source "$HOME/miniforge3/bin/activate" "cell2cell"

python3 "$work_dir/Python/202_cci_cell2cell.py" \
    --output_dir ${output_dir} \
    --meta $meta \
    --annot $annot \
    --interactions $resource \
    --nperm $nperm \
    --input_dir $sample_dir \
    --sample_id $sample_id
