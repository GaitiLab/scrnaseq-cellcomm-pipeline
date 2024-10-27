
include { RRA; COMBINE_SAMPLES } from "../nf-modules/consensus.nf"
include { AGGREGATION_CCI } from '../nf-subworkflows/aggregation_cci.nf'

workflow CCI_CONSENSUS {
    take:
    matched_cci
    metadata_rds

    main:
    RRA(
        matched_cci,
        alpha           = params.alpha,
        n_perm          = params.n_perm
    )

    COMBINE_SAMPLES(
        RRA.out.mvoted_interactions.collect(),
        RRA.out.signif_interactions.collect(),
        RRA.out.interactions_agg_rank.collect(),
        metadata            = metadata_rds,
        condition_var       = params.condition_var,
        sample_var          = params.sample_var,
        patient_var         = params.patient_var
    )

    AGGREGATION_CCI(
        mvoted_interactions     = COMBINE_SAMPLES.out.mvoted_interactions,
        interactions_agg_rank   = COMBINE_SAMPLES.out.interactions_agg_rank
    )

    emit:
    aggregation_integration = AGGREGATION_CCI.out.aggregation_integration
}
