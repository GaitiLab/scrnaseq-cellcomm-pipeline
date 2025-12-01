process SEURAT_SUBSET_OBJECT {
    label 'mem_32G'
    label 'time_30m'

    input:
    path input_file
    val sample_var
    path samplesheet

    output:
    path "${input_file.simpleName}_subset.rds", emit: rds
    path "versions.yml", emit: versions

    script:
    """
    03_subset_object.R \
    --sample_var ${sample_var} \
    --samplesheet ${samplesheet} \
    --input_file ${input_file} \
    --output_dir \${PWD} \
    --nf-process-id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 000_data
    touch "${input_file.simpleName}_subset.rds"
    touch "versions.yml"

    """
}
