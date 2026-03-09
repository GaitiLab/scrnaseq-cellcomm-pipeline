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
    ch_mvoted = params.consensus_dir
        ? channel.fromPath("${params.consensus_dir}/*interactions_mvoted.rds")
        : channel.empty()

    ch_ranked = params.consensus_dir
        ? channel.fromPath("${params.consensus_dir}/*interactions_agg_rank.rds")
        : channel.empty()

    ch_metadata = params.metadata_rds ? channel.fromPath(params.metadata_rds) : metadata_rds

    ch_combined_ranked_samples = channel.empty()
    ch_combined_mvoted_samples = channel.empty()

    if (!params.consensus_dir) {
        TAKE_CONSENSUS_ACROSS_TOOLS(matched_cci, alpha)
        TAKE_CONSENSUS_ACROSS_TOOLS.out.rds.collect().set { ch_mvoted }
        ch_versions = ch_versions.mix(TAKE_CONSENSUS_ACROSS_TOOLS.out.versions)

        RRA(matched_cci, n_perm)
        RRA.out.rds.collect().set { ch_ranked }
        ch_versions = ch_versions.mix(RRA.out.versions)
    }

    if (params.combine_samples) {
        COMBINE_RANKED_SAMPLES(ch_ranked, ch_metadata, condition_var, sample_var, patient_var, "interactions_agg_rank")
        ch_versions = ch_versions.mix(COMBINE_RANKED_SAMPLES.out.versions)
        ch_combined_ranked_samples = COMBINE_RANKED_SAMPLES.out.rds

        COMBINE_MVOTED_SAMPLES(ch_mvoted, ch_metadata, condition_var, sample_var, patient_var, "interactions_mvoted")
        ch_versions = ch_versions.mix(COMBINE_MVOTED_SAMPLES.out.versions)
        ch_combined_mvoted_samples = COMBINE_MVOTED_SAMPLES.out.rds
    }

    emit:
    ranked_rds = ch_combined_ranked_samples
    mvoted_rds = ch_combined_mvoted_samples
    versions   = ch_versions
}
