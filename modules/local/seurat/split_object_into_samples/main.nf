process SEURAT_SPLIT_OBJECT_INTO_SAMPLES {
    label 'mem_32G'
    label 'time_30m'

    input:
    path input_file
    val sample_var

    output:
    path "*.rds", emit: rds
    path "versions.yml", emit: versions

    script:
    """
    002_split_seurat_object.R \
    --input_file "${input_file}" \
    --output_dir "\${PWD}" \
    --sample_var ${sample_var} \
    --task_id ${task.process}


    """

    stub:
    """
    #!/usr/bin/env bash
    touch "split_by_${sample_var}/Sample_1.rds"
    touch "split_by_${sample_var}/Sample_2.rds"
    touch "split_by_${sample_var}/Sample_3.rds"
    touch "split_by_${sample_var}/Sample_4.rds"
    touch "split_by_${sample_var}/Sample_5.rds"
    touch "split_by_${sample_var}/Sample_6.rds"
    touch "versions.yml"

    """
}
