process RUN_CPDB {
    label 'cpdb_env'
    label 'mem_16G'
    label 'time_1h'

    input:
    tuple val(sample_id), path(barcodes), path(genes), path(matrix)
    path meta
    path interactions_db
    val annot
    val n_perm
    val min_pct

    output:
    tuple val(sample_id), path("statistical_analysis_interaction_scores__${sample_id}.txt"), path("statistical_analysis_pvalues__${sample_id}.txt"), path("statistical_analysis_significant_means__${sample_id}.txt"), path("statistical_analysis_means__${sample_id}.txt"), path("statistical_analysis_deconvoluted__${sample_id}.txt"), path("statistical_analysis_deconvoluted_percents__${sample_id}.txt"), emit: txt
    tuple val(sample_id), path("${sample_id}_counts.h5ad"), emit: h5ad
    tuple val(sample_id), path("${sample_id}_metadata.tsv"), emit: tsv

    script:
    """
    #!/usr/bin/env bash

    mkdir -p ${sample_id}

    mv ${barcodes} ${sample_id}/barcodes.tsv
    mv ${genes} ${sample_id}/genes.tsv
    mv ${matrix} ${sample_id}/matrix.mtx
    203_cci_cpdb.py \
        --input_dir \${PWD}/${sample_id} \
        --n_perm ${n_perm} \
        --interactions_db ${interactions_db} \
        --annot ${annot} \
        --sample_id ${sample_id} \
        --meta ${meta} \
        --min_pct ${min_pct} \
        --n_cores ${task.cpus} \
        --output_dir \$PWD
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "statistical_analysis_interaction_scores__${sample_id}.txt"
    touch "statistical_analysis_pvalues__${sample_id}.txt"
    touch "statistical_analysis_significant_means__${sample_id}.txt"
    touch "statistical_analysis_means__${sample_id}.txt"
    touch "statistical_analysis_deconvoluted__${sample_id}.txt
    touch "statistical_analysis_deconvoluted_percents__${sample_id}.txt"
    touch "${sample_id}_counts.h5ad"
    touch "${sample_id}_metadata.tsv"
    """
}
