process COMBINE_SAMPLES {
    label "mem_16G"
    label "time_30m"

    input:
    path "*.rds"
    path metadata
    val condition_var
    val sample_var
    val patient_var
    val suffix

    output:
    path "${suffix}.rds", emit: rds
    path "versions.yml", emit: versions

    script:
    """
    401_combine_samples.R \
    --output_dir \$PWD \
    --input_dir \$PWD \
    --metadata ${metadata} \
    --condition_var ${condition_var} \
    --sample_var ${sample_var} \
    --patient_var ${patient_var} \
    --task_id ${task.process} \
    --suffix ${suffix}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "${suffix}.rds"
    touch "versions.yml"

    """
}
