process RUN_LIANA {
    label 'mem_8G'
    label 'time_30m'

    input:
    tuple val(sample_id), path(input_file)
    path interactions_db
    val annot
    val n_perm
    val min_cells
    val min_pct

    output:
    tuple val(sample_id), path("liana__${sample_id}.rds"), emit: rds

    script:
    """
    201_cci_liana.R \
        --n_perm ${n_perm} \
        --interactions_db ${interactions_db} \
        --annot ${annot} \
        --gene_expr ${input_file} \
        --min_cells ${min_cells} \
        --min_pct ${min_pct} \
        --output_dir \${PWD}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "liana__${sample_id}.rds"
    """
}
