
process FILTER_BY_DETECTION_IN_MULTI_SAMPLES {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path input_file
    val condition_var
    val min_patients

    output:
    path "402a_filtering_detect_in_multi_samples.rds", emit: rds

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/bin/402a_filter_by_detection_in_multi_samples.R" \
    --output_dir \$PWD \
    --input_file \$PWD/${input_file} \
    --condition_var ${condition_var} \
    --min_patients ${min_patients}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "402a_filtering_detect_in_multi_samples.rds"
    """
}
