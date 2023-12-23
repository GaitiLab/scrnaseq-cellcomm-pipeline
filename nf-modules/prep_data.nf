process ADAPT_ANNOTATION {
    label 'mem64'
    label 'time_30m'

    publishDir "${projectDir}/output/${params.output_run_name}/", mode: "copy"

    input:
    path input_file
    path samples_oi

    output:
    path "000_data/seurat_annot_adapted.rds"

    when: params.do_annot

    script:
    def time_out_limit = (task.time).toSeconds() - 30

    """
    #!/usr/bin/env bash
    timeout ${time_out_limit} Rscript "${projectDir}/scripts/000a_adapt_annot.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/000_data" \
    --samples_oi "\$PWD/${samples_oi}"
    """
}

process GET_METADATA {
    label 'mem32'
    label 'time_10m'
    publishDir "${projectDir}/output/${params.output_run_name}/", mode: "copy"

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

    timeout ${time_out_limit} Rscript "${projectDir}/scripts/000b_get_metadata.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/000_data" \
    """
}

process REDUCE_SEURAT_OBJECT_SIZE {
    label 'mem32'
    label 'time_15m'

    publishDir "${projectDir}/output/${params.output_run_name}/", mode: "copy"

    input:
    path input_file

    output:
    path "000_data/${input_file.simpleName}_reduced_size.rds"

    when: (!(params.skip_reduction && file(params.sample_dir).isDirectory())) && (input_file.exists())

    script:
    def time_out_limit = (task.time).toSeconds() - 30
    """
    #!/usr/bin/env bash
    timeout ${time_out_limit} Rscript "${projectDir}/scripts/001_reduce_seurat_object_size.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/000_data"
    """
}

process SPLIT_SEURAT_OBJECT {
    label 'mem20'
    label 'time_15m'

    publishDir "${projectDir}/output/${params.output_run_name}/", mode: "copy"

    input:
    path input_file
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
    --split_varname ${split_varname}

    """

}

process PREPROCESSING {
    label 'mem4'
    label 'time_15m'

    publishDir "${projectDir}/output/${params.output_run_name}/", mode: "copy"

    input:
    path input_file
    path interactions_db
    val annot
    val min_cells
    val min_cell_types

    output:
    path "100_preprocessing/seurat/*.rds", emit: seurat_obj, optional: true
    path "100_preprocessing/mtx/*", emit: mtx_dir, optional: true

    when: (!params.skip_preprocessing) && (input_file.isFile()) && (interactions_db.isFile())

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