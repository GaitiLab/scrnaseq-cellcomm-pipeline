process FILTER_BY_DETECTION_IN_MULTI_SAMPLES {
    label "mem_16G"
    label "time_30m"

    input:
    path interactions_mvoted
    val min_patients

    output:
    path "filtering_detect_in_multi_samples.rds", emit: rds
    path "versions.yml", emit: versions

    script:
    """
    43_filter_by_detection_in_multi_samples.R \
    --output_dir \$PWD \
    --input_file ${interactions_mvoted} \
    --min_patients ${min_patients} \
    --nf-process-id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "filtering_detect_in_multi_samples.rds"
    touch "versions.yml"

    """
}
