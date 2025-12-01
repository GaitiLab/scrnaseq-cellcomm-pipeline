process CPDB_FORMAT {
    label 'time_10m'
    label 'mem_4G'

    input:
    tuple val(meta), path(interaction_scores), path(pvalues), path(significant_means), path(means), path(deconvoluted), path(deconvoluted_percents), path(ref_db)

    output:
    tuple val(meta), path("cpdb__${meta.sample_id}__postproc.rds"), emit: rds
    path "versions.yml", emit: versions

    script:
    """
    33_postproc_cpdb.R \
    --output_dir "\$PWD" \
    --sample_id ${meta.sample_id} \
    --interaction_scores ${interaction_scores} \
    --pval ${pvalues} \
    --sign_means ${significant_means} \
    --means ${means} \
    --ref_db ${ref_db} \
    --nf-process-id ${task.process}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "cpdb__${meta.sample_id}__postproc.rds"
    touch "versions.yml"

    """
}
