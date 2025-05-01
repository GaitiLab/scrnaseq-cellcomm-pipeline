process SEURAT_EXTRACT_SAMPLE {
    label 'mem_32G'
    label 'time_30m'

    input:
    tuple val(meta), path(input_file)
    val sample_var

    output:
    tuple val(meta), path("${meta.sample_id}.rds"), emit: rds
    path "versions.yml", emit: versions

    script:
    """
    003_extract_sample.R \
        --input_file ${input_file} \
        --output_dir "\${PWD}" \
        --sample_id ${meta.sample_id} \
        --task_id ${task.process} \
        --sample_var ${sample_var}

    """

    stub:
    seurat_dummy = "${meta.sample_id}.rds"
    """
    mkdir -p \$PWD/seurat
    touch "versions.yml"

    touch ${seurat_dummy}
    """
}
