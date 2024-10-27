
process REDUCE_SEURAT_OBJECT_SIZE {
    label 'mem_32G'
    label 'time_30m'
    // Commented, cause we probably do not need this object afterwards + takes a lot of space
    // publishDir params.output_dir, mode: "symlink"

    input:
    path input_file

    output:
    path "${input_file.simpleName}_reduced_size.rds", emit:rds

    script:
    """
    #!/usr/bin/env bash
    Rscript "${projectDir}/bin/001_reduce_seurat_object_size.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/000_data"
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 000_data
    touch "${input_file.simpleName}_reduced_size.rds"
    """
}
