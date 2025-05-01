process SEURAT_PREPROCESS_SAMPLE {
    label 'mem_32G'
    label 'time_30m'

    input:
    tuple val(meta), path(input_file)
    val annot
    val min_cells

    output:
    tuple val(meta), path("seurat/${meta.sample_id}.rds"), emit: rds
    tuple val(meta), path("mtx/${meta.sample_id}/barcodes.tsv"), path("mtx/${meta.sample_id}/genes.tsv"), emit: tsv
    tuple val(meta), path("mtx/${meta.sample_id}/matrix.mtx"), emit: mtx
    path "versions.yml", emit: versions

    script:
    """
    100_preprocessing.R \
        --input_file ${input_file} \
        --annot "${annot}" \
        --min_cells ${min_cells} \
        --output_dir "\$PWD" \
        --sample_id ${meta.sample_id} \
        --task_id ${task.process}

    """

    stub:
    seurat_dummy = "seurat/${meta.sample_id}.rds"
    barcodes_dummy = "mtx/${meta.sample_id}/barcodes.tsv"
    genes_dummy = "mtx/${meta.sample_id}/genes.tsv"
    matrix_dummy = "mtx/${meta.sample_id}/matrix.mtx"
    """
    mkdir -p \$PWD/mtx/${meta.sample_id}
    mkdir -p \$PWD/seurat

    touch ${seurat_dummy}
    touch ${barcodes_dummy}
    touch ${genes_dummy}
    touch ${matrix_dummy}
    touch "versions.yml"

    """
}
