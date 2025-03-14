process COMBINE_SAMPLES {
    label "mem_4G"
    label "time_10m"

    input:
    path "*__interactions_mvoted.rds"
    path "*__signif_interactions.rds"
    path "*__interactions_agg_rank.rds"
    path metadata
    val condition_var
    val sample_var
    val patient_var

    output:
    tuple path("401_samples_interactions_mvoted.rds"), path("401_samples_sign_interactions.rds"), path("401_samples_interactions_agg_rank.rds"), emit: rds
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
    --task_id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash

    mkdir -p 401_combine_samples
    touch "401_samples_interactions_mvoted.rds"
    touch "401_samples_sign_interactions.rds"
    touch "401_samples_interactions_agg_rank.rds"
    touch "versions.yml"

    """
}
