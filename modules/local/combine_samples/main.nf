process COMBINE_SAMPLES {
    label "mem_16G"
    label "time_30m"

    input:
    path paths
    path metadata
    val condition_var
    val sample_var
    val patient_var
    val suffix

    output:
    path "samples_${suffix}.rds", emit: rds
    path "versions.yml", emit: versions

    script:
    """
    42_add_metadata.R \
    --output_dir \$PWD \
    --input_dir \$PWD \
    --meta_df ${metadata} \
    --condition_var ${condition_var} \
    --sample_var ${sample_var} \
    --patient_var ${patient_var} \
    --nf-process-id ${task.process} \
    --suffix ${suffix}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "${suffix}.rds"
    touch "versions.yml"

    """
}
