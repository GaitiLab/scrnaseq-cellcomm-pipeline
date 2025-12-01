process FILTER_AGGREGATED_RESULTS {
    label "mem_16G"
    label "time_30m"

    input:
    path interactions_mvoted
    path interactions_ranked

    output:
    path "filtering_aggregated_res.rds", emit: rds
    path "versions.yml", emit: versions

    script:
    """

    44_filtering_aggregated_res.R \
    --output_dir \$PWD \
    --interactions_mvoted \$PWD/${interactions_mvoted} \
    --interactions_ranked \$PWD/${interactions_ranked} \
    --nf-process-id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "filtering_aggregated_res.rds"
    touch "versions.yml"

    """
}
