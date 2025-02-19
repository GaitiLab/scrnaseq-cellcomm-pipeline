process CREATE_SAMPLE_SHEET {
    label 'mem_16G'
    label 'time_10m'

    input:
    path input_file
    val sample_var
    val annot
    val min_cells
    val is_confident

    output:
    path "sample_sheet.csv", emit: csv

    script:
    """
    create_sample_sheet.R \
        --sample_var ${sample_var} \
        --input_file ${input_file} \
        --annot "${annot}" \
        --min_cells ${min_cells} \
        --is_confident ${is_confident} \
        --output_dir "\$PWD"
    """

    stub:
    """
    echo Sample,celltype1,celltype2 >> sample_sheet.csv


    """
}
