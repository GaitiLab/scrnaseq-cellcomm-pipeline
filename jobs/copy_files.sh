#!/usr/bin/env bash
#SBATCH -J copy_files
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

run_name="CellClass_L4_min3_types"
# project_dir="${HOME}/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions"
project_dir="/cluster/projects/gaitigroup/Users/Joan/h4h-cell-cell-interactions"

# tmp="/Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/output/CellClass_L4_min3_types/output_files.txt"
# output_dir="${project_dir}/output/${run_name}"
# while read line; do
#     parentdir=$(basename $(dirname ${line}))
#     filename=$(basename ${line})
#     # echo $parentdir
#     # echo ${parentdir}
#     if [[ $line != $output_dir/$parentdir ]]; then
#         # rm $output_dir/$parentdir/$filename
#         cp $line $output_dir/$parentdir
#     fi
# done < "$output_dir/output_files.txt"
# cp SOURCE DEST



output_dir="${project_dir}/output/${run_name}"
mkdir -p $output_dir/100_preprocessing
cd $output_dir/100_preprocessing
while read line; do
    parentdir=$(basename $(dirname ${line}))
    mkdir -p $parentdir
    cd $parentdir
    filename=$(basename ${line})
    # echo $parentdir
    # echo ${parentdir}
    # if [[ $line != $output_dir/$parentdir ]]; then
    #     # rm $output_dir/$parentdir/$filename
    cp -R $line .
    # fi
    cd ../
done < "$output_dir/preprocessing.txt"