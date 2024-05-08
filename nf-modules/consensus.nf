process CONSENSUS {
    label "mem_4G"
    label "time_30m"

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(cellchat_obj), path(liana_obj), path(cell2cell_obj), path(cpdb_obj)
    val alpha

    output:
    path("400_consensus/${sample_id}__interactions_mvoted.rds"
    ), emit: mvoted_interactions
    path("400_consensus/${sample_id}__signif_interactions.rds"
    ), emit: signif_interactions
    path("400_consensus/${sample_id}__interactions_agg_rank.rds"), emit: interactions_agg_rank

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/400_consensus.R" \
    --output_dir "\$PWD/400_consensus" \
    --sample_id ${sample_id} \
    --alpha ${alpha} \
    --cellchat_obj \$PWD/${cellchat_obj} \
    --liana_obj \$PWD/${liana_obj} \
    --cell2cell_obj \$PWD/${cell2cell_obj} \
    --cpdb_obj \$PWD/${cpdb_obj}
    """
}

process COMBINE_SAMPLES {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path "*__interactions_mvoted.rds"
    path "*__signif_interactions.rds"
    path "*__interactions_agg_rank.rds"
    path metadata
    path meta_vars_oi
    val condition_varname
    val sample_varname
    val patient_varname

    output:
    path "401_combine_samples/401_samples_interactions_mvoted.rds", emit: mvoted_interactions
    path "401_combine_samples/401_samples_sign_interactions.rds", emit: signif_interactions
    path "401_combine_samples/401_samples_interactions_agg_rank.rds", emit: interactions_agg_rank

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/401_combine_samples.R" \
    --output_dir \$PWD/401_combine_samples \
    --input_dir \$PWD \
    --metadata \$PWD/${metadata} \
    --meta_vars_oi \$PWD/${meta_vars_oi} \
    --condition_varname ${condition_varname} \
    --sample_varname ${sample_varname} \
    --patient_varname ${patient_varname}
    """
}

process AGGREGATION_BINARIZED {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path input_file
    val annot
    val condition_varname
    val min_patients

    output:
    path "402_aggregation/402a_aggregation_binarized.rds"

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/402a_aggregation_binarized.R" \
    --output_dir \$PWD/402_aggregation \
    --input_file \$PWD/${input_file} \
    --annot ${annot} \
    --condition_varname ${condition_varname} \
    --min_patients ${min_patients}
    """
}

process AGGREGATION_CONTINUOUS {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path input_file
    val condition_varname

    output:
    path "402_aggregation/402b_aggregation_continuous.rds"

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/402b_aggregation_continuous.R" \
    --output_dir \$PWD/402_aggregation \
    --input_file \$PWD/${input_file} \
    --condition_varname ${condition_varname}
    """
}

process AGGREGATION_INTEGRATION {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path interactions_agg_binarized
    path interactions_agg_continuous
    val condition_varname

    output:
    path "402_aggregation/402c_aggregation_integration.rds", emit: aggregation_integration

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/402c_aggregation_integration.R" \
    --output_dir \$PWD/402_aggregation \
    --interactions_agg_binarized \$PWD/${interactions_agg_binarized} \
    --interactions_agg_continuous \$PWD/${interactions_agg_continuous} \
    --condition_varname ${condition_varname}
    """
}


process COLLECT_RESULTS {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path interactions_agg_integration
    val condition_varname
    val alpha
    val output_name

    output:
    path "${output_name}.xlsx"

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/403_collect_results.R" \
    --output_dir \$PWD \
    --output_name ${output_name} \
    --interactions_agg_integration \$PWD/${interactions_agg_integration} \
    --condition_varname ${condition_varname} \
    --alpha ${alpha}
    """
}
