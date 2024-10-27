process formatting_cellchat {
    label 'time_10m'
    label 'mem_4G'


    input:
    tuple val(sample_id), path(input_interactions)
    path ref_db

    output:
    tuple val(sample_id), path("cellchat__${sample_id}__postproc.rds"), emit: rds

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/bin/300_postproc_cellchat.R" \
    --output_dir "\$PWD/" \
    --input_interactions \$PWD/$input_interactions \
    --ref_db \$PWD/${ref_db} \
    --sample_id ${sample_id}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cellchat__${sample_id}__postproc.rds"
    """
}
