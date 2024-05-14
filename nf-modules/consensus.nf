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

    stub: 
    """
    #!/usr/bin/env bash

    mkdir -p 400_consensus
    touch "400_consensus/${sample_id}__interactions_mvoted.rds"
    touch "400_consensus/${sample_id}__signif_interactions.rds"
    touch "400_consensus/${sample_id}__interactions_agg_rank.rds"
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
    --condition_varname ${condition_varname} \
    --sample_varname ${sample_varname} \
    --patient_varname ${patient_varname}
    """
    
    stub:
    """
    #!/usr/bin/env bash

    mkdir -p 401_combine_samples
    touch "401_combine_samples/401_samples_interactions_mvoted.rds"
    touch "401_combine_samples/401_samples_sign_interactions.rds"
    touch "401_combine_samples/401_samples_interactions_agg_rank.rds"
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

    stub: 
    """
    #!/usr/bin/env bash
    mkdir -p 402_aggregation
    touch "402_aggregation/402a_aggregation_binarized.rds"
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

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 402_aggregation
    touch "402_aggregation/402b_aggregation_continuous.rds"
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

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 402_aggregation
    touch "402_aggregation/402c_aggregation_integration.rds"
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

    stub:
    """
    #!/usr/bin/env bash
    touch "${output_name}.xlsx"
    """
}
