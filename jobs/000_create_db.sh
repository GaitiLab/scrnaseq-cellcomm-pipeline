#!/usr/bin/env bash
#SBATCH -J 200_cci_cellchat
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

echo "Activating conda environment..."
source "$HOME/miniforge3/bin/activate" "standard_env"

base_dir="${HOME}/Desktop/gaitigroup/Users"
# base_dir="/cluster/projects/gaitigroup/Users"
work_dir="$base_dir/Joan/scrnaseq-cellcomm"

output_dir="${work_dir}/001_data/interactions_db_v2"
cpdb_dir="${work_dir}/000_misc/references/cellphonedb_v5.0.0"
ref_dir="${work_dir}/000_misc/references"

Rscript "${work_dir}/scripts/011_base_liana_db.R" \
    --output_dir ${output_dir} 


Rscript "${work_dir}/scripts/012_get_cellchat_db.R" \
    --output_dir "${output_dir}"

Rscript "${work_dir}/scripts/013_get_cpdb.R" \
    --output_dir "${output_dir}" \
    --cpdb ${cpdb_dir}

Rscript "${work_dir}/scripts/014_update_liana.R" \
    --output_dir ${output_dir}

Rscript "${work_dir}/scripts/015_liana_filtering.R" \
    --output_dir ${output_dir} \
    --ref_dir ${ref_dir} \
    --interactions_dir ${output_dir}

Rscript "${work_dir}/scripts/016_update_cellchat.R" \
    --output_dir ${output_dir}

liana_db="${output_dir}/liana_db.rds"

Rscript "${work_dir}/scripts/017_update_cellphonedb.R" \
    --output_dir ${output_dir}/cellphonedb_custom \
    --cpdb ${cpdb_dir} --ref_dir "${ref_dir}" \
    --liana_db ${liana_db}


custom_cpdb_dir="${output_dir}/cellphonedb_custom"
python "${work_dir}/Python/017_update_cellphonedb.py" \
    --input_dir ${custom_cpdb_dir} 


