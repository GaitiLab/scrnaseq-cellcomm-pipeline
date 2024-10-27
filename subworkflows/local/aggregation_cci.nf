include { FILTER_BY_DETECTION_IN_MULTI_SAMPLES; AGGREGATE_SAMPLES; filtering_aggregated_results } from "../nf-modules/consensus.nf"

workflow AGGREGATION_CCI {
    take:
    mvoted_interactions
    interactions_agg_rank

    main:
    FILTER_BY_DETECTION_IN_MULTI_SAMPLES(
        input_file                      = mvoted_interactions,
        condition_var                   = params.condition_var,
        min_patients                    = params.min_patients
    )
    AGGREGATE_SAMPLES(
        input_file                      = interactions_agg_rank,
        condition_var                   = params.condition_var
    )
    filtering_aggregated_results(
        interactions_agg_binarized      = FILTER_BY_DETECTION_IN_MULTI_SAMPLES.out,
        interactions_agg_continuous     = AGGREGATE_SAMPLES.out,
        condition_var                   = params.condition_var
    )

    emit:
    aggregation_integration = filtering_aggregated_results.out.aggregation_integration

}
