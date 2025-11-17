process SEURAT_PREPROCESS_SAMPLE {
    label 'mem_32G'
    label 'time_30m'

    input:
    tuple val(meta), path(input_file)
    val annot
    val min_cells

    output:
    tuple val(meta), path("seurat/${meta.sample_id}.rds"), emit: rds
    path "versions.yml", emit: versions

    script:
    """
    10_preprocess_sample.R \
        --input_file ${input_file} \
        --annot "${annot}" \
        --min_cells ${min_cells} \
        --output_dir "\$PWD" \
        --sample_id ${meta.sample_id} \
        --nf_process_id ${task.process}

    """

    stub:
    seurat_dummy = "seurat/${meta.sample_id}.rds"

    """
    mkdir -p \$PWD/mtx/${meta.sample_id}
    mkdir -p \$PWD/seurat

    touch ${seurat_dummy}
    touch "versions.yml"

    """
}
