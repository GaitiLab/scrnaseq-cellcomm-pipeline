process GET_METADATA {
    label 'mem32'
    label 'time_10m'
    publishDir "${projectDir}/output/${params.run_name}/", mode: "copy"

    input:
    path input_file

    output:
    path "000_data/${input_file.simpleName}__metadata.rds", emit: metadata_rds
    path "000_data/${input_file.simpleName}__metadata.csv", emit: metadata_csv

    when: (!(file(params.metadata_csv).isFile() || file(params.metadata_rds).isFile())) && (input_file.isFile())

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/000_get_metadata.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/000_data" \
    """
}

process REDUCE_SEURAT_OBJECT_SIZE {
    label 'mem32'
    label 'time_15m'
    // Commented, cause we probably do not need this object afterwards + takes a lot of space
    // publishDir "${projectDir}/output/${params.run_name}/", mode: "copy"

    input:
    path input_file
    val annot

    output:
    path "000_data/${input_file.simpleName}_reduced_size.rds"

    when: (!(params.skip_reduction && file(params.sample_dir).isDirectory())) && (input_file.exists())

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash
    timeout ${time_out_limit} Rscript "${projectDir}/scripts/001_reduce_seurat_object_size.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/000_data" \
    --annot ${annot}
    """
}

process SPLIT_SEURAT_OBJECT {
    label 'mem20'
    label 'time_15m'

    publishDir "${projectDir}/output/${params.run_name}/", mode: "symlink"

    input:
    path input_file
    path downsampling_sheet
    val split_varname

    output:
    path "000_data/split_by_${split_varname}/*.rds"

    when: (!file(params.sample_dir).isDirectory()) && (input_file.isFile())

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/002_split_seurat_object.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/000_data/split_by_${split_varname}" \
    --split_varname ${split_varname} \
    --downsampling_sheet ${downsampling_sheet}

    """

}

process PREPROCESSING {
    label 'mem4'
    label 'time_15m'

    publishDir "${projectDir}/output/${params.run_name}/", mode: "copy"

    input:
    path input_file
    path interactions_db
    val annot
    val min_cells
    val min_cell_types

    output:
    path "100_preprocessing/seurat/*.rds", emit: seurat_obj, optional: true
    path "100_preprocessing/mtx/*", emit: mtx_dir, optional: true

    when: (!params.skip_preprocessing) && (input_file.isFile())

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/100_preprocessing.R" \
    --input_file \$PWD/${input_file} \
    --output_dir "\$PWD/100_preprocessing" \
    --annot "${annot}" \
    --min_cells ${min_cells} \
    --interactions_db \$PWD/${interactions_db} \
    --min_cell_types ${min_cell_types} \
    """
}