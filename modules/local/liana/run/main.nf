process LIANA_RUN {
    label 'mem_8G'
    label 'time_30m'

    input:
    tuple val(meta), path(input_file), path(interactions_db)
    val annot
    val n_perm
    val min_cells
    val min_pct

    output:
    tuple val(meta), path("liana__${meta.sample_id}.rds"), emit: rds
    path "versions.yml", emit: versions

    script:
    """
    24_cci_liana.R \
        --n_perm ${n_perm} \
        --interactions_db ${interactions_db} \
        --annot ${annot} \
        --gene_expr ${input_file} \
        --min_cells ${min_cells} \
        --min_pct ${min_pct} \
        --output_dir \${PWD} \
        --nf_process_id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "liana__${meta.sample_id}.rds"
    touch "versions.yml"

    """
}
