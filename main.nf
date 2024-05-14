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

def get_sample_id = {
    def filename_without_ext = it.simpleName
    def chunks = filename_without_ext.split( '__' )
    sample_id=chunks[1]
        return tuple(sample_id, it)
}

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
    output_dir              = file(params.output_dir)

     // Initialize empty channels
    preprocessing_mtx_dir       = Channel.empty()
    preprocessing_seurat_obj    = Channel.empty()
    matched_cci                 = Channel.empty()

    // Overwrite parameters depending on params.init_step
    if (params.init_step > 1) {
        Channel.fromPath("${output_dir}/000_data/${input_file.simpleName}__metadata.csv").set { metadata_csv }
        Channel.fromPath("${output_dir}/000_data/${input_file.simpleName}__metadata.rds").set { metadata_rds }
    }
    if (params.init_step == 2) {
        Channel.fromFilePairs("${output_dir}/100_preprocessing/mtx/*/{barcodes.tsv,genes.tsv,matrix.mtx}", type: 'file', size: -1 ) { file -> file.getParent().getName() }
            .map { id, files -> tuple(id, files).flatten() }
            .set { preprocessing_mtx_dir }
        Channel.fromPath("${output_dir}/100_preprocessing/seurat/*.rds", type : "file")
            .map { file -> tuple(file.simpleName, file)}
            .set { preprocessing_seurat_obj }
    }
    if (params.init_step == 4) {
        cci_cellchat =  Channel.fromPath("${output_dir}/300_postproc_cellchat/cellchat__*postproc.rds", type: 'file')
        | map ( get_sample_id ) 
        cci_liana = Channel.fromPath("${output_dir}/301_postproc_liana/liana__*postproc.rds", type: 'file')
        | map ( get_sample_id ) 
        cci_cell2cell = Channel.fromPath("${output_dir}/302_postproc_cell2cell/cell2cell__*.rds", type: 'file')
        | map ( get_sample_id ) 
        cci_cpdb = Channel.fromPath("${output_dir}/303_postproc_cpdb/cpdb__*postproc.rds", type: 'file')
        | map ( get_sample_id )
        matched_cci = cci_cellchat.join(cci_liana).join(cci_cell2cell).join(cci_cpdb)
    }
    cci_aggregation_integration = params.init_step == 5 ? file("${output_dir}/402_aggregation/402c_aggregation_integration.rds") : Channel.empty()
    

    println """\
    PIPELINE CONFIGURATION:
    ---- GENERAL ----------------------------------------------------------------------------------
    Project Directory   ${projectDir}
    Home Directory      ${HOME}
    
    ---- SKIP OR EXECUTE STEPS --------------------------------------------------------------------
    Skip reduction:     ${params.skip_reduction}
    Confident annot:    ${params.is_confident}
    ---- INPUTS -----------------------------------------------------------------------------------
    Input file:         ${input_file}
    Metadata csv:       ${metadata_csv}
    Metadata rds:       ${metadata_rds}

    ---- OUTPUTS ----------------------------------------------------------------------------------
    Output dir:         ${output_dir}

    ---- PRE-PROCESSING --------------------------------------------------------------------------
    Annotation:         ${params.annot}
    Min. cells:         ${params.min_cells}
    Split varname:      ${params.split_varname}

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

    // // Setup modules
    if (params.init_step == 1) {
        PREP_DATA(
            input_file                      = input_file,
        )
        PREP_DATA.out.metadata_csv.set  { metadata_csv }
        PREP_DATA.out.metadata_rds.set  { metadata_rds }
        PREP_DATA.out.mtx_dir.set       { preprocessing_mtx_dir }
        PREP_DATA.out.seurat_obj.set    { preprocessing_seurat_obj }
    }
    // preprocessing_mtx_dir.view()
    // preprocessing_seurat_obj.view()
    if (params.init_step <= 3) {
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
        matched_cci = CCI_SIMPLE_PIPELINE.out
    }
    if (params.init_step <= 4) {
        CCI_CONSENSUS(
                matched_cci                 = matched_cci,
                metadata_rds                = metadata_rds,
        )
        cci_aggregation_integration         = CCI_CONSENSUS.out.aggregation_integration
    } 
    if (params.init_step <= 5) {
        COLLECT_RESULTS(
            interactions_agg_integration    = cci_aggregation_integration, 
            condition_varname               = params.condition_varname, 
            alpha                           = params.alpha, 
            output_name                     = params.interactions_excel_name
        )
    }
}

// workflow.onComplete {

//     println ( workflow.success ? """
//         Pipeline execution summary
//         ---------------------------
//         Completed at: ${workflow.complete}
//         Duration    : ${workflow.duration}
//         Success     : ${workflow.success}
//         workDir     : ${workflow.workDir}
//         exit status : ${workflow.exitStatus}
//         """ : """
//         Failed: ${workflow.errorReport}
//         exit status : ${workflow.exitStatus}
//         """
//     )

// }