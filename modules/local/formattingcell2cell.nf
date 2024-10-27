process FORMATTING_CELL2CELL {
    label 'time_10m'
    label 'mem_4G'


    input:
    tuple val(sample_id), path(input_interactions_scores), path(input_interactions_pval)
    path ref_db

    output:
    tuple val(sample_id), path("cell2cell__${sample_id}__postproc.rds"), emit: rds

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/bin/302_postproc_cell2cell.R" \
    --output_dir "\$PWD/" \
    --input_interactions_scores \$PWD/$input_interactions_scores \
    --input_interactions_pval \$PWD/$input_interactions_pval \
    --sample_id ${sample_id} \
    --ref_db \$PWD/${ref_db}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cell2cell__${sample_id}__postproc.rds"
    """
}

