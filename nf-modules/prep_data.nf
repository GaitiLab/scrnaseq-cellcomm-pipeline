process GET_METADATA {
    label 'mem_32G'
    label 'time_10m'
    publishDir params.output_dir, mode: "copy"

    input:
    path input_file

    output:
    path "000_data/${input_file.simpleName}__metadata.rds", emit: metadata_rds
    path "000_data/${input_file.simpleName}__metadata.csv", emit: metadata_csv

    script:
    """
    #!/usr/bin/env bash
    Rscript "${projectDir}/scripts/000_get_metadata.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/000_data" \
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 000_data
    touch "000_data/${input_file.simpleName}__metadata.rds"
    touch "000_data/${input_file.simpleName}__metadata.csv"
    """   
}

process REDUCE_SEURAT_OBJECT_SIZE {
    label 'mem_32G'
    label 'time_30m'
    // Commented, cause we probably do not need this object afterwards + takes a lot of space
    publishDir params.output_dir, mode: "copy"

    input:
    path input_file
    val annot

    output:
    path "000_data/${input_file.simpleName}_reduced_size.rds"

    script:
    """
    #!/usr/bin/env bash
    Rscript "${projectDir}/scripts/001_reduce_seurat_object_size.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/000_data" \
    --annot ${annot}
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 000_data
    touch "000_data/${input_file.simpleName}_reduced_size.rds"
    """   
}

process SPLIT_SEURAT_OBJECT {
    label 'mem_32G'
    label 'time_30m'

    publishDir params.output_dir, mode: "copy"

    input:
    path input_file
    val split_varname

    output:
    path "000_data/split_by_${split_varname}/*.rds"

    script:
    """
    #!/usr/bin/env bash
    Rscript "${projectDir}/scripts/002_split_seurat_object.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/000_data/split_by_${split_varname}" \
    --split_varname ${split_varname}

    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 000_data/split_by_${split_varname}
    touch "000_data/split_by_${split_varname}/Sample_1.rds"
    touch "000_data/split_by_${split_varname}/Sample_2.rds"
    touch "000_data/split_by_${split_varname}/Sample_3.rds"
    touch "000_data/split_by_${split_varname}/Sample_4.rds"
    touch "000_data/split_by_${split_varname}/Sample_5.rds"
    touch "000_data/split_by_${split_varname}/Sample_6.rds"
    """   
}

process PREPROCESSING {
    label 'mem_8G'
    label 'time_10m'

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(input_file)
    val annot
    val min_cells
    val is_confident

    output:
    tuple val(sample_id), path("100_preprocessing/seurat/${sample_id}.rds"), emit: seurat_obj
    tuple val(sample_id), path("100_preprocessing/mtx/${sample_id}/barcodes.tsv"), 
    path("100_preprocessing/mtx/${sample_id}/genes.tsv"), path("100_preprocessing/mtx/${sample_id}/matrix.mtx"), emit: mtx_dir

    script:
    seurat_dummy="\$PWD/100_preprocessing/seurat/${sample_id}.rds"
    barcodes_dummy="\$PWD/100_preprocessing/mtx/${sample_id}/barcodes.tsv"
    genes_dummy="\$PWD/100_preprocessing/mtx/${sample_id}/genes.tsv"
    matrix_dummy="\$PWD/100_preprocessing/mtx/${sample_id}/matrix.mtx"
    """
    #!/usr/bin/env bash
    Rscript "${projectDir}/scripts/100_preprocessing.R" \
        --input_file \$PWD/${input_file} \
        --annot "${annot}" \
        --min_cells ${min_cells} \
        --is_confident ${is_confident} \
        --output_dir "\$PWD/100_preprocessing" 

    // Create empty files
    if [ -f "${seurat_dummy}" ]; then
        echo "${seurat_dummy} exists."
    else 
        echo "${seurat_dummy} does not exist."
        echo ${seurat_dummy}
        echo ${barcodes_dummy}
        echo ${genes_dummy}
        echo ${matrix_dummy}
        mkdir -p \$PWD/100_preprocessing/seurat
        mkdir -p \$PWD/100_preprocessing/mtx/${sample_id}

        touch ${seurat_dummy}
        touch ${barcodes_dummy}
        touch ${genes_dummy}
        touch ${matrix_dummy} 
    fi 
    """

    stub:
    seurat_dummy="\$PWD/100_preprocessing/seurat/${sample_id}.rds"
    barcodes_dummy="\$PWD/100_preprocessing/mtx/${sample_id}/barcodes.tsv"
    genes_dummy="\$PWD/100_preprocessing/mtx/${sample_id}/genes.tsv"
    matrix_dummy="\$PWD/100_preprocessing/mtx/${sample_id}/matrix.mtx"
    """
    mkdir -p \$PWD/100_preprocessing/mtx/${sample_id}
    mkdir -p \$PWD/100_preprocessing/seurat

    touch ${seurat_dummy}
    touch ${barcodes_dummy}
    touch ${genes_dummy}
    touch ${matrix_dummy} 
    """
}