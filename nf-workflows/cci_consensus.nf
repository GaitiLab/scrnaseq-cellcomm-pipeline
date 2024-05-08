
include { CONSENSUS; COMBINE_SAMPLES } from "../nf-modules/consensus.nf"
include { AGGREGATION_CCI } from '../nf-subworkflows/aggregation_cci.nf'

workflow CCI_CONSENSUS {
    take:
    matched_cci
    metadata_rds
    meta_vars_oi

    main: 
    CONSENSUS(
        matched_cci, 
        alpha           = params.alpha
    )

    COMBINE_SAMPLES(
        CONSENSUS.out.mvoted_interactions.collect(), 
        CONSENSUS.out.signif_interactions.collect(),
        CONSENSUS.out.interactions_agg_rank.collect(), 
        metadata            = metadata_rds, 
        meta_vars_oi        = meta_vars_oi, 
        condition_varname   = params.condition_varname, 
        sample_varname      = params.split_varname, 
        patient_varname     = params.patient_varname
    )

    AGGREGATION_CCI(
        mvoted_interactions     = COMBINE_SAMPLES.out.mvoted_interactions,
        interactions_agg_rank   = COMBINE_SAMPLES.out.interactions_agg_rank
    )

    emit: 
    aggregation_integration = AGGREGATION_CCI.out.aggregation_integration
}
