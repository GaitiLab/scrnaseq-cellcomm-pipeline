#!/usr/bin/env bash
#SBATCH -J launch_pipeline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Joan.Kant@uhn.ca
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

module load java/18

# Symlink plugins
mkdir -p $PWD/.nextflow
rm -rf $PWD/.nextflow/plugins
ln -s /cluster/home/t119972uhn/.nextflow/plugins $PWD/.nextflow/

export NXF_OFFLINE='true'

# TODO change to path to directory with the cloned Github repo.  
project_dir="/cluster/projects/gaitigroup/Users/Joan/gaitilab-scrnaseqcellcomm"

echo "PIPELINE CONFIGURATION..."
# ----  NEXTFLOW CONFIGURATION ---- #

# Work directory - all executed tasks (processes) are stored here
work_dir="${PWD}/nf-work"
# Output directory for: trace, report + timeline by NextFlow

nf_profile="apptainer,h4h,slurm"

config_file="${project_dir}/gaitilab.config"

# Create directories
mkdir -p ${work_dir}

# Test data
nextflow run ${project_dir} \
    -profile ${nf_profile} \
    -w ${work_dir} \
    -params-file "params.yml" \
    -c "${config_file}" \
    --outdir "output" -resume

echo "Done!"


