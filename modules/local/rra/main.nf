process RRA {
    label "mem_4G"
    label "time_30m"

    input:
    tuple val(meta), path(cell2cell_obj), path(cellchat_obj), path(cpdb_obj), path(liana_obj)
    val alpha
    val n_perm

    output:
    tuple path("${meta.sample_id}__interactions_mvoted.rds"), path("${meta.sample_id}__signif_interactions.rds"), path("${meta.sample_id}__interactions_agg_rank.rds"), emit: rds
    path "versions.yml", emit: versions

    script:
    """
    400_consensus_and_RRA.R \
    --output_dir \$PWD \
    --sample_id ${meta.sample_id} \
    --alpha ${alpha} \
    --cellchat_obj ${cellchat_obj} \
    --liana_obj ${liana_obj} \
    --cell2cell_obj ${cell2cell_obj} \
    --cpdb_obj ${cpdb_obj} \
    --n_perm ${n_perm} \
    --task_id ${task.process}

    """

    stub:
    """
    #!/usr/bin/env bash

    touch "${meta.sample_id}__interactions_mvoted.rds"
    touch "${meta.sample_id}__signif_interactions.rds"
    touch "${meta.sample_id}__interactions_agg_rank.rds"
    touch "versions.yml"

    """
}
