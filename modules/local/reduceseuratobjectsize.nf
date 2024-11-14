process REDUCE_SEURAT_OBJECT_SIZE {
    label 'mem_32G'
    label 'time_30m'

    input:
    path input_file

    output:
    path "${input_file.simpleName}_reduced_size.rds", emit: rds

    script:
    """
    001_reduce_seurat_object_size.R \
    --input_file ${input_file} \
    --output_dir \${PWD}
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 000_data
    touch "${input_file.simpleName}_reduced_size.rds"
    """
}
