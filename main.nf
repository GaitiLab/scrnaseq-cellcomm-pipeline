#! /usr/bin/env nextflow

/*
Workflow for cell-cell-interactions
*/

nextflow.enable.dsl=2

// Include modules
include { COLLECT_RESULTS } from './nf-modules/consensus.nf'

// Include subworkflows
include { PREP_DATA } from './nf-subworkflows/prep_data.nf'

// Include workflows
include { CCI_SIMPLE_PIPELINE } from './nf-workflows/cci.nf'
include { CCI_CONSENSUS } from './nf-workflows/cci_consensus.nf'

workflow {
    // Convert string paths
    input_file              = file(params.input_file)
    metadata_csv            = file(params.metadata_csv)
    metadata_rds            = file(params.metadata_rds)
    cellphone_db            = file(params.cellphone_db)
    cellchat_db             = file(params.cellchat_db)
    liana_db                = file(params.liana_db)
    cell2cell_db            = file(params.cell2cell_db)
    ref_db                  = file(params.ref_db)
    meta_vars_oi            = file(params.meta_vars_oi)
    output_dir              = file(params.output_dir)

    // Create empty channels (initalization)
    preprocessing_mtx_dir = Channel.fromPath("${output_dir}/100_preprocessing/mtx", type = "dir")
    preprocessing_seurat_obj = Channel.fromPath("${output_dir}/100_preprocessing/seurat/*.rds", type = "file")

    println """\
    PIPELINE CONFIGURATION:
    ---- GENERAL ----------------------------------------------------------------------------------

    ---- SKIP OR EXECUTE STEPS --------------------------------------------------------------------
    Skip reduction:     ${params.skip_reduction}

    ---- INPUTS -----------------------------------------------------------------------------------
    Input file:         ${input_file}
    Metadata csv:       ${metadata_csv}
    Metadata rds:       ${metadata_rds}
    Meta vars oi:       ${meta_vars_oi}

    ---- OUTPUTS ----------------------------------------------------------------------------------
    Output dir:         ${output_dir}

    ---- PRE-PROCESSING --------------------------------------------------------------------------
    Annotation:         ${params.annot}
    Min. cells:         ${params.min_cells}
    Split varname:      ${params.split_varname}
    Min. cell types:    ${params.min_cell_types}

    ---- DATABASES --------------------------------------------------------------------------------
    CellphoneDB:        ${cellphone_db}
    CellChatDB:         ${cellchat_db}
    LianaDB:            ${liana_db}
    Cell2cell:          ${cell2cell_db}
    ReferenceDB:        ${ref_db}
    
    ---- CCI CONFIG / FORMATTING ------------------------------------------------------------------
    Annotation:         ${params.annot}
    Min. pct:           ${params.min_pct}
    Alpha (sign level): ${params.alpha}
    N permutations:     ${params.n_perm}
    Condition:          ${params.condition_varname}
    Min.patients:       ${params.min_patients}
    Patient varname:    ${params.patient_varname}
    """.stripIndent()

    // Setup modules
    if (params.init_step == 1) {
            PREP_DATA(
                input_file                  = input_file,
                liana_db                    = liana_db
            )
            PREP_DATA.out.metadata_csv.set  { metadata_csv }
            PREP_DATA.out.metadata_rds.set  { metadata_rds }
            PREP_DATA.out.mtx_dir.set       { preprocessing_mtx_dir }
            PREP_DATA.out.seurat_obj.set    { preprocessing_seurat_obj }
    }

    CCI_SIMPLE_PIPELINE( 
        preprocessing_seurat_obj        = preprocessing_seurat_obj,
        preprocessing_mtx_dir           = preprocessing_mtx_dir,
        metadata_csv                    = metadata_csv,
        cellchat_db                     = cellchat_db,
        cell2cell_db                    = cell2cell_db,
        liana_db                        = liana_db,
        cellphone_db                    = cellphone_db,
        ref_db                          = ref_db,
        output_dir                      = output_dir
    )

    CCI_CONSENSUS(
        matched_cci                     = CCI_SIMPLE_PIPELINE.out.matched_cci,
        metadata_rds                    = metadata_rds,
        meta_vars_oi                    = meta_vars_oi
    )

    COLLECT_RESULTS(
        interactions_agg_integration    = CCI_CONSENSUS.out.aggregation_integration, 
        condition_varname               = params.condition_varname, 
        alpha                           = params.alpha, 
        output_name                     = params.interactions_excel_name
    )
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )

}