process RRA {
    label "mem_4G"
    label "time_30m"


    input:
    tuple val(sample_id), path(cellchat_obj), path(liana_obj), path(cell2cell_obj), path(cpdb_obj)
    val alpha
    val n_perm

    output:
    tuple
    path("${sample_id}__interactions_mvoted.rds"
    ),
    path("${sample_id}__signif_interactions.rds"
    ),
    path("${sample_id}__interactions_agg_rank.rds"), emit: rds

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/bin/400_consensus_and_RRA.R" \
    --output_dir \$PWD \
    --sample_id ${sample_id} \
    --alpha ${alpha} \
    --cellchat_obj \$PWD/${cellchat_obj} \
    --liana_obj \$PWD/${liana_obj} \
    --cell2cell_obj \$PWD/${cell2cell_obj} \
    --cpdb_obj \$PWD/${cpdb_obj} \
    --n_perm ${n_perm}
    """

    stub:
    """
    #!/usr/bin/env bash

    touch "${sample_id}__interactions_mvoted.rds"
    touch "${sample_id}__signif_interactions.rds"
    touch "${sample_id}__interactions_agg_rank.rds"
    """
}
