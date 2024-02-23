#!/usr/bin/env bash
# Copy configuration files nf-config folder
source_dir="nf-config-local"
target_dir="nf-config"

# Create target directory if not existing
mkdir -p ${target_dir}

# Using gcp on MacOS
for current_path in $source_dir/*; do 
    # echo $current_path
    abs_path=$(readlink $current_path)
    filename=$(basename $current_path)

    if [[ $abs_path == "" ]]
        then
            # echo "No symlink found for ${filename}..."
            gcp -f $current_path "nf-config"
            continue
    else 
        # echo "Symlink found for ${filename}..."
        # echo "True path ${abs_path}..."
        # echo $abs_path
        gcp -f --preserve --remove-destination "${abs_path}" "${target_dir}/${filename}"
    fi 
done 
