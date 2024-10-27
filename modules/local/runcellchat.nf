process RUN_CELLCHAT {
    label 'mem_32G'
    label 'time_12h'

    input:
    tuple val(sample_id), path(input_file)
    path interactions_db
    val annot
    val n_perm
    val min_cells

    output:
    tuple val(sample_id), path("cellchat__${sample_id}.rds"), path "cellchat__${sample_id}__raw_obj.rds", emit:rds

    script:
    """
    #!/usr/bin/env bash
    Rscript "${projectDir}/bin/200_cci_cellchat.R" \
        --n_perm ${n_perm} \
        --interactions_db \$PWD/${interactions_db} \
        --annot ${annot} \
        --gene_expr \$PWD/${input_file} \
        --min_cells ${min_cells} \
        --output_dir "\$PWD" \
        --n_cores ${task.cpus}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cellchat__${sample_id}.rds"
    touch "cellchat__${sample_id}__raw_obj.rds"
    """
}
