process CELLCHAT_RUN {
    label 'mem_64G'
    label 'time_12h'

    input:
    tuple val(meta), path(input_file), path(interactions_db)
    val annot
    val n_perm
    val min_cells

    output:
    tuple val(meta), path("cellchat__${meta.sample_id}__raw_obj.rds"), emit: rds
    path "versions.yml", emit: versions

    script:
    """
    22_cci_cellchat.R \
        --n_perm ${n_perm} \
        --interactions_db_path ${interactions_db} \
        --annot ${annot} \
        --gene_expr_path ${input_file} \
        --min_cells ${min_cells} \
        --output_dir "\$PWD" \
        --n_cores ${task.cpus} \
        --nf-process-id ${task.process} \
        --sample_id ${meta.sample_id}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cellchat__${meta.sample_id}__raw_obj.rds"
    touch "versions.yml"

    """
}
