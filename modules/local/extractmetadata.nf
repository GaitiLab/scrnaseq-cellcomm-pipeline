process EXTRACT_METADATA {
    label 'mem_32G'
    label 'time_10m'

    input:
    path input_file

    output:
    path "${input_file.simpleName}__metadata.rds", emit: rds
    path "${input_file.simpleName}__metadata.csv", emit: csv

    script:
    """
    000_get_metadata.R \
    --input_file "${input_file}" \
    --output_dir "\$PWD" \
    """

    stub:
    """
    #!/usr/bin/env bash

    touch "${input_file.simpleName}__metadata.rds"
    touch "${input_file.simpleName}__metadata.csv"
    """
}
