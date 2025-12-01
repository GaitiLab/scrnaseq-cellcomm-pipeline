include { RRA                                       } from "../../../modules/local/rra"
include { COMBINE_SAMPLES as COMBINE_MVOTED_SAMPLES ; COMBINE_SAMPLES as COMBINE_RANKED_SAMPLES } from '../../../modules/local/combine_samples'
include { TAKE_CONSENSUS_ACROSS_TOOLS               } from '../../../modules/local/take_consensus_across_tools'

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
    ch_versions = channel.empty()
    ch_mvoted = params.interactions
        ? channel.fromPath("${params.interactions}/*interactions_mvoted.rds")
        : channel.empty()

    ch_ranked = params.interactions
        ? channel.fromPath("${params.interactions}/*interactions_agg_rank.rds")
        : channel.empty()

    ch_metadata = params.metadata_rds ? channel.fromPath(params.metadata_rds) : metadata_rds

    if (!params.interactions) {
        TAKE_CONSENSUS_ACROSS_TOOLS(matched_cci, alpha)
        TAKE_CONSENSUS_ACROSS_TOOLS.out.rds.collect().set { ch_mvoted }
        ch_versions = ch_versions.mix(TAKE_CONSENSUS_ACROSS_TOOLS.out.versions)

        RRA(matched_cci, n_perm)
        RRA.out.rds.collect().set { ch_ranked }
        ch_versions = ch_versions.mix(RRA.out.versions)
    }

    COMBINE_RANKED_SAMPLES(ch_ranked, ch_metadata, condition_var, sample_var, patient_var, "interactions_agg_rank")
    ch_versions = ch_versions.mix(COMBINE_RANKED_SAMPLES.out.versions)

    COMBINE_MVOTED_SAMPLES(ch_mvoted, ch_metadata, condition_var, sample_var, patient_var, "interactions_mvoted")
    ch_versions = ch_versions.mix(COMBINE_MVOTED_SAMPLES.out.versions)

    emit:
    ranked_rds = COMBINE_RANKED_SAMPLES.out.rds
    mvoted_rds = COMBINE_MVOTED_SAMPLES.out.rds
    versions   = ch_versions
}
