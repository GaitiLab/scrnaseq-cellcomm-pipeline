process FORMATTING_CPDB {
    label 'time_10m'
    label 'mem_4G'

    input:
    tuple val(sample_id), path(interaction_scores),
    path(pvalues), path(significant_means), path(means)
    path ref_db

    output:
    tuple val(sample_id), path("cpdb__${sample_id}__postproc.rds"), emit:rds

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/bin/303_postproc_cellphonedb.R" \
    --output_dir "\$PWD" \
    --sample_id ${sample_id} \
    --interaction_scores \$PWD/${interaction_scores} \
    --pval \$PWD/${pvalues} \
    --sign_means \$PWD/${significant_means} \
    --means \$PWD/${means} \
    --ref_db \$PWD/${ref_db}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cpdb__${sample_id}__postproc.rds"
    """
}
