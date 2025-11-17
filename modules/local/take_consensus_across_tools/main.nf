process TAKE_CONSENSUS_ACROSS_TOOLS {
    label "mem_4G"
    label "time_30m"

    input:
    tuple val(meta), path(cell2cell_obj), path(cellchat_obj), path(cpdb_obj), path(liana_obj)
    val alpha

    output:
    tuple val(meta), path("${meta.sample_id}__interactions_mvoted.rds"), emit: rds
    path "versions.yml", emit: versions

    script:
    """
    41_take_consensus_across_tools.R \
    --output_dir \$PWD \
    --sample_id ${meta.sample_id} \
    --alpha ${alpha} \
    --cellchat_obj ${cellchat_obj} \
    --liana_obj ${liana_obj} \
    --cell2cell_obj ${cell2cell_obj} \
    --cpdb_obj ${cpdb_obj} \
    --nf_process_id ${task.process}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "${meta.sample_id}__interactions_mvoted.rds"
    touch "versions.yml"

    """
}
