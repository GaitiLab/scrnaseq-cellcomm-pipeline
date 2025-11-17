process CPDB_RUN {
    label 'cpdb_env'
    label 'mem_16G'
    label 'time_1h'

    input:
    tuple val(meta), path(barcodes), path(genes), path(matrix), path(meta_path), path(interactions_db)
    val annot
    val n_perm
    val min_pct

    output:
    tuple val(meta), path("statistical_analysis_interaction_scores__${meta.sample_id}.txt"), path("statistical_analysis_pvalues__${meta.sample_id}.txt"), path("statistical_analysis_significant_means__${meta.sample_id}.txt"), path("statistical_analysis_means__${meta.sample_id}.txt"), path("statistical_analysis_deconvoluted__${meta.sample_id}.txt"), path("statistical_analysis_deconvoluted_percents__${meta.sample_id}.txt"), emit: txt
    tuple val(meta), path("${meta.sample_id}_counts.h5ad"), emit: h5ad
    tuple val(meta), path("${meta.sample_id}_metadata.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env bash

    mkdir -p ${meta.sample_id}

    mv ${barcodes} ${meta.sample_id}/barcodes.tsv
    mv ${genes} ${meta.sample_id}/genes.tsv
    mv ${matrix} ${meta.sample_id}/matrix.mtx
    23_cci_cpdb.py \
        --input_dir \${PWD}/${meta.sample_id} \
        --n_perm ${n_perm} \
        --interactions_db ${interactions_db} \
        --annot ${annot} \
        --sample_id ${meta.sample_id} \
        --meta ${meta_path} \
        --min_pct ${min_pct} \
        --n_cores ${task.cpus} \
        --output_dir \$PWD \
        --nf_process_id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "statistical_analysis_interaction_scores__${meta.sample_id}.txt"
    touch "statistical_analysis_pvalues__${meta.sample_id}.txt"
    touch "statistical_analysis_significant_means__${meta.sample_id}.txt"
    touch "statistical_analysis_means__${meta.sample_id}.txt"
    touch "statistical_analysis_deconvoluted__${meta.sample_id}.txt
    touch "statistical_analysis_deconvoluted_percents__${meta.sample_id}.txt"
    touch "${meta.sample_id}_counts.h5ad"
    touch "${meta.sample_id}_metadata.tsv"
    touch "versions.yml"

    """
}
