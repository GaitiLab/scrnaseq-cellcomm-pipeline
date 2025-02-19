process AGGREGATE_SAMPLES {
    label "mem_4G"
    label "time_10m"

    input:
    tuple path(_interactions_mvoted), path(_signif_interactions), path(interactions_agg_rank)
    val condition_var

    output:
    path "402b_aggregation_samples.rds", emit: rds

    script:
    """
    402b_aggregation_samples.R \
    --output_dir \$PWD \
    --input_file ${interactions_agg_rank} \
    --condition_var ${condition_var}
    """

    stub:
    """
    #!/usr/bin/env bash

    touch "402b_aggregation_samples.rds"
    """
}
