process consensus_and_rra {
    label "mem_4G"
    label "time_30m"

    publishDir params.output_dir, mode: "copy"

    input:
    tuple val(sample_id), path(cellchat_obj), path(liana_obj), path(cell2cell_obj), path(cpdb_obj)
    val alpha

    output:
    path("400_consensus_and_RRA/${sample_id}__interactions_mvoted.rds"
    ), emit: mvoted_interactions
    path("400_consensus_and_RRA/${sample_id}__signif_interactions.rds"
    ), emit: signif_interactions
    path("400_consensus_and_RRA/${sample_id}__interactions_agg_rank.rds"), emit: interactions_agg_rank

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/400_consensus_and_RRA.R" \
    --output_dir "\$PWD/400_consensus_and_RRA" \
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

    mkdir -p 400_consensus_and_RRA
    touch "400_consensus_and_RRA/${sample_id}__interactions_mvoted.rds"
    touch "400_consensus_and_RRA/${sample_id}__signif_interactions.rds"
    touch "400_consensus_and_RRA/${sample_id}__interactions_agg_rank.rds"
    """
}

process combine_samples {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path "*__interactions_mvoted.rds"
    path "*__signif_interactions.rds"
    path "*__interactions_agg_rank.rds"
    path metadata
    val condition_var
    val sample_var
    val patient_var

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
    --condition_var ${condition_var} \
    --sample_var ${sample_var} \
    --patient_var ${patient_var}
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

process filtering_detect_in_multi_samples {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path input_file
    val condition_var
    val min_patients

    output:
    path "402_aggregation_and_filtering/402a_filtering_detect_in_multi_samples.rds"

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/402a_filtering_detect_in_multi_samples.R" \
    --output_dir \$PWD/402_aggregation_and_filtering \
    --input_file \$PWD/${input_file} \
    --condition_var ${condition_var} \
    --min_patients ${min_patients}
    """

    stub: 
    """
    #!/usr/bin/env bash
    mkdir -p 402_aggregation_and_filtering
    touch "402_aggregation_and_filtering/402a_filtering_detect_in_multi_samples.rds"
    """
}

process aggregation_samples {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path input_file
    val condition_var

    output:
    path "402_aggregation_and_filtering/402b_aggregation_samples.rds"

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/402b_aggregation_samples.R" \
    --output_dir \$PWD/402_aggregation_and_filtering \
    --input_file \$PWD/${input_file} \
    --condition_var ${condition_var}
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 402_aggregation_and_filtering
    touch "402_aggregation_and_filtering/402b_aggregation_samples.rds"
    """
}

process filtering_aggregated_results {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path interactions_agg_binarized
    path interactions_agg_continuous
    val condition_var

    output:
    path "402_aggregation_and_filtering/402c_filtering_aggregated_res.rds", emit: aggregation_integration

    script:
    """
    #!/usr/bin/env bash

    Rscript "${projectDir}/scripts/402c_filtering_aggregated_res.R" \
    --output_dir \$PWD/402_aggregation_and_filtering \
    --interactions_agg_binarized \$PWD/${interactions_agg_binarized} \
    --interactions_agg_continuous \$PWD/${interactions_agg_continuous} \
    --condition_var ${condition_var}
    """

    stub:
    """
    #!/usr/bin/env bash
    mkdir -p 402_aggregation_and_filtering
    touch "402_aggregation_and_filtering/402c_filtering_aggregated_res.rds"
    """
}


process save_as_xlsx {
    label "mem_4G"
    label "time_10m"

    publishDir params.output_dir, mode: "copy"

    input:
    path interactions_agg_integration
    val condition_var
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
    --condition_var ${condition_var} \
    --alpha ${alpha}
    """

    stub:
    """
    #!/usr/bin/env bash
    touch "${output_name}.xlsx"
    """
}
