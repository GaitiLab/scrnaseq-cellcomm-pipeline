process CELL2CELL_RUN {
    label 'mem_64G'
    label 'time_24h'

    input:
    tuple val(meta), path(barcodes), path(genes), path(matrix), path(meta_path), path(interactions_db)
    val annot
    val n_perm

    output:
    tuple val(meta), path("cell2cell__${meta.sample_id}.pickle"), emit: pickle
    tuple val(meta), path("cell2cell__${meta.sample_id}__pvalues.csv"), path("cell2cell__${meta.sample_id}__interaction_scores.csv"), emit: csv
    path "versions.yml", emit: versions

    script:
    """

    202_cci_cell2cell.py \
        --input_dir \$PWD \
        --n_perm ${n_perm} \
        --interactions_db ${interactions_db} \
        --annot ${annot} \
        --sample_id ${meta.sample_id} \
        --meta ${meta_path} \
        --output_dir "\$PWD" \
        --nf_process_id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cell2cell__${meta.sample_id}.pickle"
    touch "cell2cell__${meta.sample_id}__pvalues.csv"
    touch "cell2cell__${meta.sample_id}__interaction_scores.csv"
    touch "versions.yml"

    """
}
