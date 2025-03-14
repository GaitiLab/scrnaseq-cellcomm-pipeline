process CELLCHAT_RUN {
    label 'mem_32G'
    label 'time_12h'

    input:
    tuple val(meta), path(input_file), path(interactions_db)
    val annot
    val n_perm
    val min_cells

    output:
    tuple val(meta), path("cellchat__${meta.sample_id}.rds"), path("cellchat__${meta.sample_id}__raw_obj.rds"), emit: rds
    path "versions.yml", emit: versions

    script:
    """
    200_cci_cellchat.R \
        --n_perm ${n_perm} \
        --interactions_db ${interactions_db} \
        --annot ${annot} \
        --gene_expr ${input_file} \
        --min_cells ${min_cells} \
        --output_dir "\$PWD" \
        --n_cores ${task.cpus} \
        --task_id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cellchat__${meta.sample_id}.rds"
    touch "cellchat__${meta.sample_id}__raw_obj.rds"
    touch "versions.yml"

    """
}
