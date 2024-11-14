include { RRA             } from "../../../modules/local/rra.nf"
include { COMBINE_SAMPLES } from '../../../modules/local/combinesamples.nf'


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
    RRA(
        matched_cci,
        alpha,
        n_perm
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

    COMBINE_SAMPLES(result.mvoted.collect(), result.signif.collect(), result.agg_rank.collect(), metadata_rds, condition_var, sample_var, patient_var)

    emit:
    rds = COMBINE_SAMPLES.out.rds
}
