include { AGGREGATE_SAMPLES                    } from '../../../modules/local/aggregate_samples'
include { FILTER_BY_DETECTION_IN_MULTI_SAMPLES } from '../../../modules/local/filter_by_detection_in_multi_samples'
include { FILTER_AGGREGATED_RESULTS            } from '../../../modules/local/filter_aggregated_results'
include { UTILS_SAVE_AS_XLSX                   } from '../../../modules/local/utils/save_as_xlsx/main.nf'

workflow AGGREGATION {
    take:
    consensus_objects
    condition_var
    min_patients

    main:

    ch_versions = Channel.empty()

    FILTER_BY_DETECTION_IN_MULTI_SAMPLES(
        consensus_objects,
        condition_var,
        min_patients,
    )
    ch_versions = ch_versions.mix(FILTER_BY_DETECTION_IN_MULTI_SAMPLES.out.versions)

    AGGREGATE_SAMPLES(consensus_objects, condition_var)
    ch_versions = ch_versions.mix(AGGREGATE_SAMPLES.out.versions)


    FILTER_AGGREGATED_RESULTS(
        FILTER_BY_DETECTION_IN_MULTI_SAMPLES.out.rds,
        AGGREGATE_SAMPLES.out.rds,
        condition_var,
    )
    ch_versions = ch_versions.mix(FILTER_AGGREGATED_RESULTS.out.versions)
    ch_aggregated_cci = FILTER_AGGREGATED_RESULTS.out.rds

    if (params.export_as_excel) {
        UTILS_SAVE_AS_XLSX(ch_aggregated_cci, params.condition_var, params.alpha, params.interactions_excel_name)
        ch_versions = ch_versions.mix(UTILS_SAVE_AS_XLSX.out.versions)
    }

    emit:
    rds      = ch_aggregated_cci
    versions = ch_versions
}
