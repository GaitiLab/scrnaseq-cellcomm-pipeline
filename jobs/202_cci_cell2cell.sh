#!/usr/bin/env bash
#SBATCH -J 202_cci_cell2cell
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
##SBATCH --partition=himem
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=01:00:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out


nperm=10
annot="CellClass_L4"

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/scrnaseq-cellcomm
sample_id="6467_solid_core"
resource="${work_dir}/001_data/interactions_db_v2/liana_db.csv"

meta="${work_dir}/output/CellClass_L4_min3_types_rerun/000_data/seurat_annot_adapted__metadata.csv"
sample_dir="${work_dir}/output/CellClass_L4_min3_types_rerun/100_preprocessing/mtx/6467_solid_core"
output_dir="${work_dir}/final_output/"

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
EOF