#! /usr/bin/env nextflow

/*
Workflow for cell-cell-interactions

*/

nextflow.enable.dsl=2

include { ADAPT_ANNOTATION; GET_METADATA; REDUCE_SEURAT_OBJECT_SIZE; PREPROCESSING; SPLIT_SEURAT_OBJECT } from "./nf-modules/prep_data.nf"
include { INFER_CELLCHAT; INFER_LIANA; INFER_CELL2CELL; INFER_CPDB } from "./nf-modules/infer_interactions.nf"
include { POSTPROCESSING_CELLCHAT; POSTPROCESSING_LIANA; POSTPROCESSING_CELL2CELL; POSTPROCESSING_CPDB } from "./nf-modules/filtering.nf"
include { CONSENSUS; COMBINE_SAMPLES; AGGREGATION_PATIENT; AGGREGATION_SAMPLE } from "./nf-modules/consensus.nf"

workflow {
    // Convert string paths
    input_file = file(params.input_file)
    sample_dir = file(params.sample_dir)
    metadata_csv = file(params.metadata_csv)
    metadata_rds = file(params.metadata_rds)
    cellphone_db = file(params.cellphone_db)
    cellchat_db = file(params.cellchat_db)
    liana_db = file(params.liana_db)
    liana_db_csv = file(params.liana_db_csv)
    ref_db = file(params.ref_db)
    meta_vars_oi = file(params.meta_vars_oi)
    samples_oi = file(params.samples_oi)

    println """\
    PIPELINE CONFIGURATION:
    ---- GENERAL ----------------------------------------------------------------------------------
    Run name:           ${params.output_run_name}
    Approach:           ${params.approach}

    ---- SKIP OR EXECUTE STEPS --------------------------------------------------------------------
    Do annot:           ${params.do_annot}
    Skip reduction:     ${params.skip_reduction}
    Skip preprocessing: ${params.skip_preprocessing}

    ---- INPUTS -----------------------------------------------------------------------------------
    Input file:         ${input_file}
    Sample dir:         ${sample_dir}
    Metadata csv:       ${metadata_csv}
    Metadata rds:       ${metadata_rds}
    Meta vars oi:       ${meta_vars_oi}
    Samples oi:         ${samples_oi}
    
    ---- OUTPUTS ----------------------------------------------------------------------------------
    Output dir:         "${projectDir}/${params.output_run_name}"

    ---- PRE-PROCESSING --------------------------------------------------------------------------
    Annotation:         ${params.annot}
    Min. cells:         ${params.min_cells}
    First N celltypes:  ${params.first_n}
    Split varname:      ${params.split_varname}
    Samples of interest:${samples_oi}
    Min. cell types:    ${params.min_cell_types}

    ---- DATABASES --------------------------------------------------------------------------------
    CellphoneDB:        ${cellphone_db}
    CellChatDB:         ${cellchat_db}
    LianaDB:            ${liana_db}
    LianaDB (C2C):      ${liana_db_csv}
    ReferenceDB:        ${ref_db}
    
    ---- CCI CONFIG / FORMATTING ------------------------------------------------------------------
    Annotation:         ${params.annot}
    Min. pct:           ${params.min_pct}
    Alpha (sign level): ${params.alpha}
    N permutations:     ${params.n_perm}
    Condition:          ${params.condition_varname}
    Min.patients:       ${params.min_patients}
    Patient varname:    ${params.patient_varname}

    ---- AGGREGATION -----------------------------------------------------------------------------
    Aggregate samples:  ${params.aggregate_samples}
    Aggregate patients: ${params.aggregate_patients}


    """.stripIndent()



    if(params.approach >= 1) {
        // OPTIONAL: Re-annotate / Annotate 
        // - input file provided
        // - do_annot = false    
        ADAPT_ANNOTATION(input_file = input_file, samples_oi = samples_oi)
        input_file = (params.do_annot) ? ADAPT_ANNOTATION.out : input_file

        // Grabbing metadata from input file 
        // - if metadata files aren't provided
        GET_METADATA(input_file = input_file)

        metadata_csv = metadata_csv.isFile() ? metadata_csv : GET_METADATA.out.metadata_csv
        metadata_rds = metadata_rds.isFile() ? metadata_rds : GET_METADATA.out.metadata_rds

        // OPTIONAL: Reduce size of seurat object to reduce memory usage for inferring interactinos. 
        REDUCE_SEURAT_OBJECT_SIZE(input_file = input_file)
        input_file2 = (params.skip_reduction) ? input_file : REDUCE_SEURAT_OBJECT_SIZE.out

        // Split seurat object into samples
        SPLIT_SEURAT_OBJECT(input_file = input_file2, split_varname = params.split_varname)
        samples = (file(params.sample_dir).isDirectory() && !file(input_file2).isFile()) ? Channel.fromPath("${sample_dir}/*.rds", type: 'file') : SPLIT_SEURAT_OBJECT.out.flatten()
        
        PREPROCESSING(input_file = samples, interactions_db = liana_db,
        annot = params.annot,
        min_cells = params.min_cells,
        min_cell_types = params.min_cell_types)

        preprocessing_mtx_dir = params.skip_preprocessing ? Channel.fromPath("${sample_dir}/mtx/*", type: 'dir') : PREPROCESSING.out.mtx_dir.flatten()
        preprocessing_seurat_obj = params.skip_preprocessing ? Channel.fromPath("${sample_dir}/seurat/*.rds", type: 'file') : PREPROCESSING.out.seurat_obj.flatten()
    }
    // INFERRING INTERACTIONS
    if (params.approach >= 2 ) {
        INFER_CELLCHAT(preprocessing_seurat_obj, interactions_db = cellchat_db, annot = params.annot, n_perm = params.n_perm)

        INFER_LIANA(preprocessing_seurat_obj, interactions_db = liana_db, annot = params.annot, n_perm = params.n_perm)

        INFER_CELL2CELL(input_dir = preprocessing_mtx_dir,
        meta = metadata_csv, interactions_db = liana_db_csv, annot = params.annot, n_perm = params.n_perm)

        INFER_CPDB(input_dir = preprocessing_mtx_dir,
        meta = metadata_csv, interactions_db = cellphone_db, annot = params.annot, n_perm = params.n_perm, min_pct = params.min_pct, alpha = params.alpha)

    }
    // POST-PROCESSING INTERACTIONS
    if (params.approach >= 3) {
        POSTPROCESSING_CELLCHAT(INFER_CELLCHAT.out.cellchat_obj, interactions_db = ref_db)

        POSTPROCESSING_LIANA(INFER_LIANA.out.liana_obj, ref_db = ref_db)

        POSTPROCESSING_CELL2CELL(INFER_CELL2CELL.out.cell2cell_obj, ref_db = ref_db)

        POSTPROCESSING_CPDB(INFER_CPDB.out.cpdb_obj, interactions_db = ref_db)
    }
    // MERGE INTERACTIONS based on sample id
    if (params.approach >= 4) {
        combined_objects = POSTPROCESSING_CELLCHAT.out.join(POSTPROCESSING_LIANA.out, by: 0).join(POSTPROCESSING_CELL2CELL.out, by: 0).join(POSTPROCESSING_CPDB.out, by: 0)
        // // Take consensus - sample wise
        CONSENSUS(combined_objects, alpha = params.alpha)
    }
    if (params.approach >= 5) {
        COMBINE_SAMPLES(CONSENSUS.out.mvoted_interactions.collect(), CONSENSUS.out.signif_interactions.collect(), 
        metadata = metadata_rds, meta_vars_oi = meta_vars_oi, condition_varname = params.condition_varname, sample_varname = params.split_varname, patient_varname = params.patient_varname)
    }
    if (params.approach >= 6) {
        if(params.aggregate_samples) {
            AGGREGATION_SAMPLE(COMBINE_SAMPLES.out.mvoted_interactions, metadata = metadata_rds, min_cells = params.min_cells, min_frac_samples = params.min_frac_samples, annot = params.annot, condition_varname = params.condition_varname, sample_varname = params.split_varname )
        }
        if(params.aggregate_patients) {
            AGGREGATION_PATIENT(COMBINE_SAMPLES.out.mvoted_interactions, annot = params.annot, condition_varname = params.condition_varname, min_patients = params.min_patients)
        }
    } 
}