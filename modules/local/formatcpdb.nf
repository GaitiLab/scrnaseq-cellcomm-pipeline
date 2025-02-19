process FORMAT_CPDB {
    label 'time_10m'
    label 'mem_4G'

    input:
    tuple val(sample_id), path(interaction_scores), path(pvalues), path(significant_means), path(means), path(deconvoluted), path(deconvoluted_percents)
    path ref_db

    output:
    tuple val(sample_id), path("cpdb__${sample_id}__postproc.rds"), emit: rds

    script:
    """
    303_postproc_cellphonedb.R \
    --output_dir "\$PWD" \
    --sample_id ${sample_id} \
    --interaction_scores ${interaction_scores} \
    --pval ${pvalues} \
    --sign_means ${significant_means} \
    --means ${means} \
    --ref_db ${ref_db}

    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cpdb__${sample_id}__postproc.rds"
    """
}
