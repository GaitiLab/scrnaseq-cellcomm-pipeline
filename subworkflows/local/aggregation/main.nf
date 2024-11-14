include { AGGREGATE_SAMPLES                    } from '../../../modules/local/aggregatesamples.nf'
include { FILTER_BY_DETECTION_IN_MULTI_SAMPLES } from '../../../modules/local/filterbydetectioninmultisamples.nf'
include { FILTER_AGGREGATED_RESULTS            } from '../../../modules/local/filteraggregatedresults.nf'

workflow AGGREGATION {
    take:
    consensus_objects
    condition_var    
    min_patients     

    main:

    FILTER_BY_DETECTION_IN_MULTI_SAMPLES(
        consensus_objects,
        condition_var,
        min_patients
    )

    AGGREGATE_SAMPLES(consensus_objects, condition_var)

    FILTER_AGGREGATED_RESULTS(
        FILTER_BY_DETECTION_IN_MULTI_SAMPLES.out.rds,
        AGGREGATE_SAMPLES.out.rds,
        condition_var
    )

    emit:
    rds = FILTER_AGGREGATED_RESULTS.out.rds
}
