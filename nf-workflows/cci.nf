include { INFER_CCI } from "../nf-subworkflows/infer_cci.nf"
include { POSTPROCESSING_CCI } from "../nf-subworkflows/postprocessing_cci.nf"

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

    main:
    INFER_CCI( 
        preprocessing_seurat_obj    = preprocessing_seurat_obj,
        preprocessing_mtx_dir       = preprocessing_mtx_dir,
        metadata_csv                = metadata_csv,
        cellchat_db                 = cellchat_db,
        cell2cell_db                = cell2cell_db,
        liana_db                    = liana_db,
        cellphone_db                = cellphone_db
    )

    POSTPROCESSING_CCI( 
        liana_obj                   = INFER_CCI.out.cci_liana,
        cell2cell_obj               = INFER_CCI.out.cci_cell2cell,
        cpdb_obj                    = INFER_CCI.out.cci_cpdb, 
        cellchat_obj                = INFER_CCI.out.cci_cellchat,
        ref_db                      = ref_db
    )

    emit: 
    matched_cci                     = POSTPROCESSING_CCI.out.cci_cellchat
                                        .join(
                                        POSTPROCESSING_CCI.out.cci_liana, by: 0)
                                        .join(
                                        POSTPROCESSING_CCI.out.cci_cell2cell, by: 0)
                                        .join(
                                        POSTPROCESSING_CCI.out.cci_cpdb, by: 0)
}