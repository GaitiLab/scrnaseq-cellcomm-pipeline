process UTILS_COMBINE_SAMPLES_METADATA {
    label 'mem_32G'
    label 'time_10m'

    input:
    path input_dir

    output:
    path "metadata.rds", emit: rds
    path "metadata.csv", emit: csv
    path "versions.yml", emit: versions

    script:
    """
    05_combineSamplesMetadata.R \
    --input_dir "\$PWD" \
    --output_dir "\$PWD" \
    --nf-process-id ${task.process}
    """

    stub:
    """
    #!/usr/bin/env bash

    touch "metadata.rds"
    touch "metadata.csv"
    touch "versions.yml"

    """
}
