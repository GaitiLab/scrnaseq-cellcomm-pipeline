process UTILS_SAVE_AS_XLSX {
    label "mem_4G"
    label "time_10m"

    input:
    path interactions_agg_integration
    val alpha
    val output_name

    output:
    path "${output_name}.xlsx"
    path "versions.yml", emit: versions

    script:
    """
    45_collect_results.R \
    --output_dir \$PWD \
    --output_name ${output_name} \
    --interactions_agg_integration ${interactions_agg_integration} \
    --alpha ${alpha} \
    --nf-process-id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "${output_name}.xlsx"
    touch "versions.yml"
    """
}
