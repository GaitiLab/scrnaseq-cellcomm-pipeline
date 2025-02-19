#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/scrnaseqcellcomm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/scrnaseqcellcomm
    Website: https://nf-co.re/scrnaseqcellcomm
    Slack  : https://nfcore.slack.com/channels/scrnaseqcellcomm
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SCRNASEQCELLCOMM } from './workflows/scrnaseqcellcomm/main.nf'

//
// WORKFLOW: Run main nf-core/scrnaseqcellcomm analysis pipeline
//
workflow NFCORE_SCRNASEQCELLCOMM {

    // Convert string paths
    input_file              = file(params.input_file)
    metadata_csv            = file(params.metadata_csv)
    metadata_rds            = file(params.metadata_rds)
    cellphone_db            = file(params.cellphone_db)
    cellchat_db             = file(params.cellchat_db)
    liana_db                = file(params.liana_db)
    cell2cell_db            = file(params.cell2cell_db)
    ref_db                  = file(params.ref_db)

    // // Initialize empty channels
    // preprocessing_mtx_dir       = Channel.empty()
    // preprocessing_seurat_obj    = Channel.empty()
    // matched_cci                 = Channel.empty()

    SCRNASEQCELLCOMM (
        input_file,
        metadata_csv,
        metadata_rds,
        params.annot,
        params.min_cells,
        params.sample_var,
        cellphone_db,
        cellchat_db,
        liana_db,
        cell2cell_db,
        ref_db,
        params.min_pct,
        params.n_perm,
        params.condition_var,
        params.patient_var,
        params.is_confident,
        params.min_patients,
        params.alpha,
        params.interactions_excel_name,
        params.skip_reduction
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_SCRNASEQCELLCOMM ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
