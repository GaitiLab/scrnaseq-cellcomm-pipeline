process SAMPLE_PREPROCESSING {
    label 'mem_8G'
    label 'time_10m'

    input:
    tuple val(sample_id), path(input_file)
    val annot
    val min_cells
    val is_confident

    output:
    tuple val(sample_id), path("seurat/${sample_id}.rds"), emit: rds
    tuple val(sample_id), path("mtx/${sample_id}/barcodes.tsv"), path("mtx/${sample_id}/genes.tsv"), emit: tsv
    tuple val(sample_id), path("mtx/${sample_id}/matrix.mtx"), emit: mtx

    script:
    seurat_dummy = "\$PWD/seurat/${sample_id}.rds"
    barcodes_dummy = "\$PWD/mtx/${sample_id}/barcodes.tsv"
    genes_dummy = "\$PWD/mtx/${sample_id}/genes.tsv"
    matrix_dummy = "\$PWD/mtx/${sample_id}/matrix.mtx"
    """
    100_preprocessing.R \
        --input_file ${input_file} \
        --annot "${annot}" \
        --min_cells ${min_cells} \
        --is_confident ${is_confident} \
        --output_dir "\$PWD" \
        --sample_id ${sample_id}
    """

    stub:
    seurat_dummy = "seurat/${sample_id}.rds"
    barcodes_dummy = "mtx/${sample_id}/barcodes.tsv"
    genes_dummy = "mtx/${sample_id}/genes.tsv"
    matrix_dummy = "mtx/${sample_id}/matrix.mtx"
    """
    mkdir -p \$PWD/mtx/${sample_id}
    mkdir -p \$PWD/seurat

    touch ${seurat_dummy}
    touch ${barcodes_dummy}
    touch ${genes_dummy}
    touch ${matrix_dummy}
    """
}
