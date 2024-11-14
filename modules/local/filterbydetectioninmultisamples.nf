process FILTER_BY_DETECTION_IN_MULTI_SAMPLES {
    label "mem_4G"
    label "time_10m"
    publishDir params.output_dir, mode: "copy"

    input:
    tuple path(interactions_mvoted), path(_signif_interactions), path(_interactions_agg_rank)
    val condition_var
    val min_patients

    output:
    path "402a_filtering_detect_in_multi_samples.rds", emit: rds

    script:
    """
    402a_filter_by_detection_in_multi_samples.R \
    --output_dir \$PWD \
    --input_file ${interactions_mvoted} \
    --condition_var ${condition_var} \
    --min_patients ${min_patients}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "402a_filtering_detect_in_multi_samples.rds"
    """
}
