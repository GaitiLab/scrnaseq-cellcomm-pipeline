process COLLECT_RESULTS_AS_XLSX {
    label "mem_4G"
    label "time_10m"

    input:
    path interactions_agg_integration
    val condition_var
    val alpha
    val output_name

    output:
    path "${output_name}.xlsx"

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/bin/403_collect_results.R" \
    --output_dir \$PWD \
    --output_name ${output_name} \
    --interactions_agg_integration \$PWD/${interactions_agg_integration} \
    --condition_var ${condition_var} \
    --alpha ${alpha}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "${output_name}.xlsx"
    """
}
