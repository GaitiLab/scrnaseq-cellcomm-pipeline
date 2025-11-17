process AGGREGATE_SAMPLES {
    label "mem_16G"
    label "time_30m"

    input:
    tuple path(_interactions_mvoted), path(interactions_agg_rank)

    output:
    path "aggregation_samples.rds.rds", emit: rds
    path "versions.yml", emit: versions

    script:
    """
    43_aggregation_samples.R \
    --output_dir \$PWD \
    --input_file ${interactions_agg_rank} \
    --nf_process_id ${task.process}
    """

    stub:
    """
    #!/usr/bin/env bash

    touch "aggregation_samples.rds.rds"
    touch "versions.yml"
    """
}
