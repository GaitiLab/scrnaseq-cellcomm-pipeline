process LIANA_FORMAT {
    label 'time_10m'
    label 'mem_4G'

    input:
    tuple val(meta), path(input_interactions), path(ref_db)

    output:
    tuple val(meta), path("liana__${meta.sample_id}__postproc.rds"), emit: rds
    path "versions.yml", emit: versions

    script:

    """

    301_postproc_liana.R \
    --output_dir "\$PWD" \
    --input_interactions ${input_interactions} \
    --sample_id ${meta.sample_id} \
    --ref_db ${ref_db} \
    --task_id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "liana__${meta.sample_id}__postproc.rds"
    touch "versions.yml"

    """
}
