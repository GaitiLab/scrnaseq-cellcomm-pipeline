include { POSTPROCESSING_CCI } from "../nf-subworkflows/postprocessing_cci.nf"
include { CELLCHAT; CPDB; LIANA; CELL2CELL } from "../nf-subworkflows/infer_cci.nf"

// Extract sample ID from path, e
def get_sample_id = {
    def filename_without_ext = it.simpleName
    def chunks = filename_without_ext.split( '__' )
    sample_id=chunks[1]
        return tuple(sample_id, it)
}

workflow CCI_SIMPLE_PIPELINE {
    take:
    preprocessing_seurat_obj   
    preprocessing_mtx_dir     
    metadata_csv                
    cellchat_db
    cell2cell_db
    liana_db
    cellphone_db
    ref_db
    output_dir

    main:
    // Initialize parameters
    matched_cci                     = Channel.empty()
    cci_cellchat                    = Channel.empty()
    cci_liana                       = Channel.empty()
    cci_cell2cell                   = Channel.empty()
    cci_cpdb                        = Channel.empty()

    // Only keep samples with at least 2 cell types
    preprocessing_mtx_dir               = preprocessing_mtx_dir | branch { it ->
                                            TRUE: it[1].size() > 0
                                            FALSE: it[1].size() <= 0
                                        }
    preprocessing_mtx_dir.TRUE.set { preprocessing_mtx_dir }

    preprocessing_seurat_obj               = preprocessing_seurat_obj | branch { it ->
                                            TRUE: it[1].size() > 0
                                            FALSE: it[1].size() <= 0
                                    }
    preprocessing_seurat_obj.TRUE.set { preprocessing_seurat_obj }
    if (params.init_step <= 2) {
        CELLCHAT(
            preprocessing_seurat_obj    = preprocessing_seurat_obj,
            cellchat_db                 = cellchat_db,
            ref_db                      = ref_db
        )

        LIANA(
            preprocessing_seurat_obj    = preprocessing_seurat_obj,
            liana_db                    = liana_db,
            ref_db                      = ref_db
        )

        CELL2CELL(
            preprocessing_mtx_dir       = preprocessing_mtx_dir,
            metadata_csv                = metadata_csv,
            cell2cell_db                = cell2cell_db,
            ref_db                      = ref_db
        )

        CPDB(
            preprocessing_mtx_dir       = preprocessing_mtx_dir,
            metadata_csv                = metadata_csv,
            cellphone_db                = cellphone_db,
            ref_db                      = ref_db
        )

        matched_cci = CELLCHAT.out.join(LIANA.out).join(CELL2CELL.out).join(CPDB.out)
    } else if (params.init_step == 3) {
        cci_cellchat    = Channel.fromPath(
            "${output_dir}/200_cci_cellchat/cellchat__*.rds", 
            type: 'file'
            )
            | map ( get_sample_id ) 
            | groupTuple( sort:true ) 
            | map { id, files -> tuple(id, files[0]) }
        cci_liana = Channel.fromPath("${output_dir}/201_cci_liana/liana__*.rds", type: 'file')
                                        | map ( get_sample_id )
        cci_cell2cell   =  Channel.fromPath(
            "${output_dir}/202_cci_cell2cell/cell2cell__*.csv", 
            type: 'file')
            | map ( get_sample_id )
        // cci_cpdb = Channel.fromPath( 
        //     [
        //         "${output_dir}/203_cci_cpdb/statistical_analysis_interaction_scores__*.txt", 
        //         "${output_dir}/203_cci_cpdb/statistical_analysis_pvalues__*.txt",
        //         "${output_dir}/203_cci_cpdb/statistical_analysis_significant_means__*.txt",
        //         "${output_dir}/203_cci_cpdb/statistical_analysis_means__*.txt"
        //     ]) 
        //     | map ( get_sample_id ) 
        //     | groupTuple 
        //     | map { id, files -> tuple(id, files[0], files[1], files[2], files[3]) } 

        cci_cpdb = Channel.fromFilePairs( "${output_dir}/203_cci_cpdb/{statistical_analysis_interaction_scores,statistical_analysis_pvalues,statistical_analysis_significant_means,statistical_analysis_means}__*.txt")
            .map { id, files -> tuple(id, files).flatten() }

        POSTPROCESSING_CCI( 
            liana_obj                   = cci_liana,
            cell2cell_obj               = cci_cell2cell,
            cpdb_obj                    = cci_cpdb, 
            cellchat_obj                = cci_cellchat,
            ref_db                      = ref_db
        )
        matched_cci = POSTPROCESSING_CCI.out.cci_cellchat
                    .join(
                    POSTPROCESSING_CCI.out.cci_liana, by: 0)
                    .join(
                    POSTPROCESSING_CCI.out.cci_cell2cell, by: 0)
                    .join(
                    POSTPROCESSING_CCI.out.cci_cpdb, by: 0)
    }
    emit: 
    matched_cci
}