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
    ch_results = params.interactions
        ? Channel.fromPath("${params.interactions}/*.rds").branch { x ->
            mvoted: x.name.endsWith("interactions_mvoted.rds")
            signif: x.name.endsWith("signif_interactions.rds")
            agg_rank: x.name.endsWith("interactions_agg_rank.rds")
        }
        : Channel.empty()
    ch_metadata = params.metadata_rds ? Channel.fromPath(params.metadata_rds) : metadata_rds


    ch_results.mvoted.view()


    if (!params.interactions) {
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
            .set { ch_results }
        ch_versions = ch_versions.mix(RRA.out.versions)
    }

    COMBINE_SAMPLES(ch_results.mvoted.collect(), ch_results.signif.collect(), ch_results.agg_rank.collect(), ch_metadata, condition_var, sample_var, patient_var)
    ch_versions = ch_versions.mix(COMBINE_SAMPLES.out.versions)

    emit:
    rds      = COMBINE_SAMPLES.out.rds
    versions = ch_versions
}
