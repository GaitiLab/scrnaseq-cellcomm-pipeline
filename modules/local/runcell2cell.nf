process RUN_CELL2CELL {
    label 'mem_64G'
    label 'time_24h'

    input:
    tuple val(sample_id), path(barcodes), path(genes), path(matrix)
    path meta
    path interactions_db
    val annot
    val n_perm

    output:
    tuple val(sample_id), path("cell2cell__${sample_id}.pickle"), emit: pickle
    tuple val(sample_id), path("cell2cell__${sample_id}__pvalues.csv"), path("cell2cell__${sample_id}__interaction_scores.csv"), emit: csv

    script:
    """

    202_cci_cell2cell.py \
        --input_dir \$PWD \
        --n_perm ${n_perm} \
        --interactions_db ${interactions_db} \
        --annot ${annot} \
        --sample_id ${sample_id} \
        --meta ${meta} \
        --output_dir "\$PWD"
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cell2cell__${sample_id}.pickle"
    touch "cell2cell__${sample_id}__pvalues.csv"
    touch "cell2cell__${sample_id}__interaction_scores.csv"
    """
}
