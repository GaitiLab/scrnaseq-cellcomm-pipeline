include { filtering_detect_in_multi_samples; aggregation_samples; filtering_aggregated_results } from "../nf-modules/consensus.nf"

workflow AGGREGATION_CCI {
    take: 
    mvoted_interactions
    interactions_agg_rank

    main:
    filtering_detect_in_multi_samples(
        input_file                      = mvoted_interactions, 
        annot                           = params.annot, 
        condition_var               = params.condition_var, 
        min_patients                    = params.min_patients
    )
    aggregation_samples(
        input_file                      = interactions_agg_rank, 
        condition_var               = params.condition_var
    )
    filtering_aggregated_results(
        interactions_agg_binarized      = filtering_detect_in_multi_samples.out, 
        interactions_agg_continuous     = aggregation_samples.out, 
        condition_var               = params.condition_var
    )

    emit: 
    aggregation_integration = filtering_aggregated_results.out.aggregation_integration

}