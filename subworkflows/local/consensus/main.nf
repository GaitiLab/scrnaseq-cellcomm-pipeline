include { RRA             } from "../../../modules/local/rra"
include { COMBINE_SAMPLES } from '../../../modules/local/combine_samples'


workflow CONSENSUS {
    take:
    matched_cci
    metadata_rds
    alpha
    n_perm
    condition_var
    sample_var
    patient_var

    main:
    ch_versions = Channel.empty()

    RRA(
        matched_cci,
        alpha,
        n_perm,
    )
    RRA.out.rds
        .collect()
        .flatten()
        .branch { x ->
            mvoted: x.name.endsWith("interactions_mvoted.rds")
            signif: x.name.endsWith("signif_interactions.rds")
            agg_rank: x.name.endsWith("interactions_agg_rank.rds")
        }
        .set { result }

    ch_versions = ch_versions.mix(RRA.out.versions)

    COMBINE_SAMPLES(result.mvoted.collect(), result.signif.collect(), result.agg_rank.collect(), metadata_rds, condition_var, sample_var, patient_var)
    ch_versions = ch_versions.mix(COMBINE_SAMPLES.out.versions)

    emit:
    rds      = COMBINE_SAMPLES.out.rds
    versions = ch_versions
}
