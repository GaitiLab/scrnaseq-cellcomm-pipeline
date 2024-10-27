

process AGGREGATE_SAMPLES {
    label "mem_4G"
    label "time_10m"

    input:
    path input_file
    val condition_var

    output:
    path "402b_aggregation_samples.rds", emit:rds

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/bin/402b_aggregation_samples.R" \
    --output_dir \$PWD \
    --input_file \$PWD/${input_file} \
    --condition_var ${condition_var}
    """

    stub:
    """
    #!/usr/bin/env bash

    touch "402b_aggregation_samples.rds"
    """
}
