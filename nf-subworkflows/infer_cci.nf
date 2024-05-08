include { INFER_CELLCHAT; INFER_LIANA; INFER_CELL2CELL; INFER_CPDB } from "../nf-modules/infer_interactions.nf"

workflow INFER_CCI {
    take: 
    preprocessing_seurat_obj
    preprocessing_mtx_dir
    metadata_csv
    cellchat_db
    cell2cell_db
    liana_db
    cellphone_db
    
    main: 
    INFER_CELLCHAT(
        preprocessing_seurat_obj, 
        interactions_db     = cellchat_db, 
        annot               = params.annot,
        n_perm              = params.n_perm
    )

    INFER_LIANA(
        preprocessing_seurat_obj, 
        interactions_db     = liana_db, 
        annot               = params.annot, 
        n_perm              = params.n_perm
    )

    INFER_CELL2CELL(
        input_dir           = preprocessing_mtx_dir,
        meta                = metadata_csv, 
        interactions_db     = cell2cell_db, 
        annot               = params.annot, 
        n_perm              = params.n_perm
    )

    INFER_CPDB(
        input_dir           = preprocessing_mtx_dir,
        meta                = metadata_csv, 
        interactions_db     = cellphone_db, 
        annot               = params.annot, 
        n_perm              = params.n_perm, 
        min_pct             = params.min_pct, 
    )

    emit: 
    cci_cellchat            = INFER_CELLCHAT.out.cellchat_obj
    cci_liana               = INFER_LIANA.out.liana_obj
    cci_cell2cell           = INFER_CELL2CELL.out.cell2cell_obj
    cci_cpdb                = INFER_CPDB.out.cpdb_obj


}