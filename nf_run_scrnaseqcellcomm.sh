#!/usr/bin/env bash
#SBATCH -J launch_pipeline-spotlight
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
##SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

export NXF_OFFLINE='true'
export NXF_DEFAULT_DSL=2

module load java/18
module load apptainer

# ${nf_exec_path} run ${PWD} -profile apptainer,h4h,debug \
#     -c "gaitilab.config" --outdir "output"

# nextflow run ${PWD} -profile apptainer,h4h -params-file "nf-params.yml" -c "gaitilab.config" -resume \
#     --outdir "output_test"


# Pipeline for testing Xenium slides
# nextflow run ${PWD} -profile apptainer,h4h,slurm -params-file "nf-params-xenium.yml" -c "gaitilab.config" -resume \
#     --outdir "output_test"


# nextflow run ${PWD} -profile apptainer,h4h,slurm -params-file "nf-params-tcga.yml" -c "gaitilab.config" -resume \
#     --outdir "output"

nextflow run ${PWD} -profile conda,slurm -params-file "nf-params.yml" --outdir "output" -c "gaitilab.config"
