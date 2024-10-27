process FORMATTING_LIANA {
    label 'time_10m'
    label 'mem_4G'


    input:
    tuple val(sample_id), path(input_interactions)
    path ref_db

    output:
    tuple val(sample_id), path("liana__${sample_id}__postproc.rds"), emit:rds

    script:

    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/bin/301_postproc_liana.R" \
    --output_dir "\$PWD" \
    --input_interactions \$PWD/$input_interactions \
    --sample_id ${sample_id} \
    --ref_db \$PWD/${ref_db}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "liana__${sample_id}__postproc.rds"
    """
}

