include { AGGREGATION_BINARIZED; AGGREGATION_CONTINUOUS; AGGREGATION_INTEGRATION } from "../nf-modules/consensus.nf"

workflow AGGREGATION_CCI {
    take: 
    mvoted_interactions
    interactions_agg_rank

    main:
    AGGREGATION_BINARIZED(
        input_file                      = mvoted_interactions, 
        annot                           = params.annot, 
        condition_varname               = params.condition_varname, 
        min_patients                    = params.min_patients
    )
    AGGREGATION_CONTINUOUS(
        input_file                      = interactions_agg_rank, 
        condition_varname               = params.condition_varname
    )
    AGGREGATION_INTEGRATION(
        interactions_agg_binarized      = AGGREGATION_BINARIZED.out, 
        interactions_agg_continuous     = AGGREGATION_CONTINUOUS.out, 
        condition_varname               = params.condition_varname
    )

    emit: 
    aggregation_integration = AGGREGATION_INTEGRATION.out.aggregation_integration

}