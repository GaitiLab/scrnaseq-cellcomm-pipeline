/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREP_DATA               } from '../../subworkflows/local/prep_data/main.nf'
include { CELLCHAT                } from "../../subworkflows/local/cellchat/main.nf"
include { CPDB                    } from "../../subworkflows/local/cpdb/main.nf"
include { LIANA                   } from '../../subworkflows/local/liana/main.nf'
include { CELL2CELL               } from '../../subworkflows/local/cell2cell/main.nf'
include { CONSENSUS               } from '../../subworkflows/local/consensus/main.nf'
include { AGGREGATION             } from '../../subworkflows/local/aggregation/main.nf'
include { COLLECT_RESULTS_AS_XLSX } from '../../modules/local/collectresultsasxlsx.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCRNASEQCELLCOMM {
    take:
    input_file             
    metadata_csv           
    metadata_rds           
    annot                  
    min_cells              
    sample_var             
    cellphone_db           
    cellchat_db            
    liana_db               
    cell2cell_db           
    ref_db                 
    min_pct                
    n_perm                 
    condition_var          
    patient_var            
    is_confident           
    min_patients           
    alpha                  
    interactions_excel_name
    skip_reduction         

    main:

    def scrnaseqcellcomm_modules = params.scrnaseqcellcomm_modules ? params.scrnaseqcellcomm_modules.split(',').collect { it.trim().toLowerCase() } : []

    if (scrnaseqcellcomm_modules.contains("prep_data")) {
        PREP_DATA(
            input_file,
            sample_var,
            annot,
            min_cells,
            is_confident,
            skip_reduction
        )

        metadata_csv = PREP_DATA.out.metadata_csv
        metadata_rds = PREP_DATA.out.metadata_rds
        mtx_dir_prepped = PREP_DATA.out.mtx_dir
        seurat_obj_prepped = PREP_DATA.out.seurat_obj
    }

    if (scrnaseqcellcomm_modules.contains("run_cci")) {
        CELL2CELL(
            mtx_dir_prepped,
            metadata_csv,
            cell2cell_db,
            ref_db,
            annot,
            n_perm
        )

        CELLCHAT(
            seurat_obj_prepped,
            cellchat_db,
            ref_db,
            annot,
            n_perm,
            min_cells
        )

        CPDB(
            mtx_dir_prepped,
            metadata_csv,
            cellphone_db,
            ref_db,
            annot,
            n_perm,
            min_pct
        )

        LIANA(
            seurat_obj_prepped,
            liana_db,
            ref_db,
            annot,
            n_perm,
            min_cells,
            min_pct
        )
        matched_cci = CELL2CELL.out.rds.join(CELLCHAT.out.rds).join(CPDB.out.rds).join(LIANA.out.rds)
    }

    if (scrnaseqcellcomm_modules.contains("consensus")) {
        CONSENSUS(
            matched_cci,
            metadata_rds,
            alpha,
            n_perm,
            condition_var,
            sample_var,
            patient_var
        )

        consensus_ch = CONSENSUS.out.rds
    }

    if (scrnaseqcellcomm_modules.contains("aggregation")) {
        AGGREGATION(
            consensus_ch,
            condition_var,
            min_patients
        )

        aggregation_ch = AGGREGATION.out.rds
    }

    if (scrnaseqcellcomm_modules.contains("export_as_excel")) {
        COLLECT_RESULTS_AS_XLSX(aggregation_ch, condition_var, alpha, interactions_excel_name)
    }
}
