
include { consensus_and_rra; combine_samples } from "../nf-modules/consensus.nf"
include { AGGREGATION_CCI } from '../nf-subworkflows/aggregation_cci.nf'

workflow CCI_CONSENSUS {
    take:
    matched_cci
    metadata_rds

    main: 
    consensus_and_rra(
        matched_cci, 
        alpha           = params.alpha,
        n_perm          = params.n_perm
    )

    combine_samples(
        consensus_and_rra.out.mvoted_interactions.collect(), 
        consensus_and_rra.out.signif_interactions.collect(),
        consensus_and_rra.out.interactions_agg_rank.collect(), 
        metadata            = metadata_rds, 
        condition_var       = params.condition_var, 
        sample_var          = params.sample_var, 
        patient_var         = params.patient_var
    )

    AGGREGATION_CCI(
        mvoted_interactions     = combine_samples.out.mvoted_interactions,
        interactions_agg_rank   = combine_samples.out.interactions_agg_rank
    )

    emit: 
    aggregation_integration = AGGREGATION_CCI.out.aggregation_integration
}
