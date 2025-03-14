process UTILS_CREATE_SAMPLESHEET {
    label 'mem_16G'
    label 'time_10m'

    input:
    path input_file
    val sample_var
    val annot
    val min_cells

    output:
    path "samplesheet.csv", emit: csv
    path "versions.yml", emit: versions

    script:
    """
    create_sample_sheet.R \
        --sample_var ${sample_var} \
        --input_file ${input_file} \
        --annot "${annot}" \
        --min_cells ${min_cells} \
        --output_dir "\$PWD" \
        --task_id ${task.process}

    """

    stub:
    """
    echo Sample,celltype1,celltype2 >> samplesheet.csv
    touch "versions.yml"


    """
}
