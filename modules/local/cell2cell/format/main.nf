process CELL2CELL_FORMAT {
    label 'time_10m'
    label 'mem_4G'

    input:
    tuple val(meta), path(input_interactions_pval), path(input_interactions_scores), path(ref_db)

    output:
    tuple val(meta), path("cell2cell__${meta.sample_id}__postproc.rds"), emit: rds
    path "versions.yml", emit: versions

    script:
    """
    31_postproc_cell2cell.R \
    --output_dir "\${PWD}" \
    --input_interactions_scores ${input_interactions_scores} \
    --input_interactions_pval ${input_interactions_pval} \
    --sample_id ${meta.sample_id} \
    --ref_db ${ref_db} \
    --nf-process-id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cell2cell__${meta.sample_id}__postproc.rds"
    touch "versions.yml"

    """
}
