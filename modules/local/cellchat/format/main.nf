process CELLCHAT_FORMAT {
    label 'time_10m'
    label 'mem_4G'

    input:
    tuple val(meta), path(input_interactions), path(ref_db)

    output:
    tuple val(meta), path("cellchat__${meta.sample_id}__postproc.rds"), emit: rds
    path "versions.yml", emit: versions

    script:
    """
    32_postproc_cellchat.R \
    --output_dir "\$PWD" \
    --input_interactions ${input_interactions} \
    --ref_db ${ref_db} \
    --sample_id ${meta.sample_id} \
    --nf-process-id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cellchat__${meta.sample_id}__postproc.rds"
    touch "versions.yml"

    """
}
