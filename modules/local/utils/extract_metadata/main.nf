process UTILS_EXTRACT_METADATA {
    label 'mem_32G'
    label 'time_10m'

    input:
    path input_file

    output:
    path "${input_file.simpleName}__metadata.rds", emit: rds
    path "${input_file.simpleName}__metadata.csv", emit: csv
    path "versions.yml", emit: versions

    script:
    """
    00_get_metadata.R \
    --input_file "${input_file}" \
    --output_dir "\$PWD" \
    --nf_process_id ${task.process}
    """

    stub:
    """
    #!/usr/bin/env bash

    touch "${input_file.simpleName}__metadata.rds"
    touch "${input_file.simpleName}__metadata.csv"
    touch "versions.yml"

    """
}
