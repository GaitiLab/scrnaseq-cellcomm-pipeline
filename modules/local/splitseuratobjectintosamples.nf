
process SPLIT_SEURAT_OBJECT_INTO_SAMPLES {
    label 'mem_32G'
    label 'time_30m'


    input:
    path input_file
    val sample_var

    output:
    path "split_by_${sample_var}/*.rds", emit:rds

    script:
    """
    #!/usr/bin/env bash
    Rscript "${projectDir}/bin/002_split_seurat_object.R" \
    --input_file "\$PWD/${input_file}" \
    --output_dir "\$PWD/split_by_${sample_var}" \
    --sample_var ${sample_var}

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
    """
}
