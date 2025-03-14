process FILTER_AGGREGATED_RESULTS {
    label "mem_4G"
    label "time_10m"

    input:
    path interactions_agg_binarized
    path interactions_agg_continuous
    val condition_var

    output:
    path "402c_filtering_aggregated_res.rds", emit: rds
    path "versions.yml", emit: versions

    script:
    """

    402c_filtering_aggregated_res.R \
    --output_dir \$PWD \
    --interactions_agg_binarized \$PWD/${interactions_agg_binarized} \
    --interactions_agg_continuous \$PWD/${interactions_agg_continuous} \
    --condition_var ${condition_var} \
    --task_id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "402c_filtering_aggregated_res.rds"
    touch "versions.yml"

    """
}
