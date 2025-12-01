process SEURAT_CONVERT_TO_MTX {
    label 'mem_32G'
    label 'time_30m'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("mtx/${meta.sample_id}/barcodes.tsv"), path("mtx/${meta.sample_id}/genes.tsv"), emit: tsv
    tuple val(meta), path("mtx/${meta.sample_id}/matrix.mtx"), emit: mtx
    path "versions.yml", emit: versions

    script:
    """
    11_convert_seurat_to_mtx.R \
        --input_file ${input_file} \
        --output_dir "\$PWD" \
        --sample_id ${meta.sample_id} \
        --nf-process-id ${task.process}

    """

    stub:
    barcodes_dummy = "mtx/${meta.sample_id}/barcodes.tsv"
    genes_dummy = "mtx/${meta.sample_id}/genes.tsv"
    matrix_dummy = "mtx/${meta.sample_id}/matrix.mtx"
    """
    mkdir -p \$PWD/mtx/${meta.sample_id}
    mkdir -p \$PWD/seurat

    touch ${barcodes_dummy}
    touch ${genes_dummy}
    touch ${matrix_dummy}
    touch "versions.yml"

    """
}
