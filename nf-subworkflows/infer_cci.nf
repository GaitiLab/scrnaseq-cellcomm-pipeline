include { INFER_CELLCHAT; INFER_LIANA; INFER_CELL2CELL; INFER_CPDB } from "../nf-modules/infer_interactions.nf"
include { POSTPROCESSING_CELLCHAT; POSTPROCESSING_LIANA; POSTPROCESSING_CELL2CELL; POSTPROCESSING_CPDB } from "../nf-modules/filtering.nf"

workflow CELLCHAT {
    take: 
    preprocessing_seurat_obj
    cellchat_db
    ref_db 

    main: 
    INFER_CELLCHAT(
        preprocessing_seurat_obj, 
        interactions_db     = cellchat_db, 
        annot               = params.annot,
        n_perm              = params.n_perm,
        min_cells           = params.min_cells
    )

    POSTPROCESSING_CELLCHAT(
        INFER_CELLCHAT.out.cellchat_obj, 
        interactions_db     = ref_db
    )

    emit: 
    POSTPROCESSING_CELLCHAT.out
}

workflow CPDB {
    take: 
    preprocessing_mtx_dir
    metadata_csv
    cellphone_db
    ref_db
    
    main: 
    INFER_CPDB(
        preprocessing_mtx_dir,
        meta                = metadata_csv, 
        interactions_db     = cellphone_db, 
        annot               = params.annot, 
        n_perm              = params.n_perm, 
        min_pct             = params.min_pct, 
    )

    POSTPROCESSING_CPDB(
        INFER_CPDB.out.cpdb_obj, 
        interactions_db     = ref_db
    )

    emit: 
    POSTPROCESSING_CPDB.out
}

workflow LIANA {
    take: 
    preprocessing_seurat_obj
    liana_db
    ref_db
    
    main: 
    INFER_LIANA(
        preprocessing_seurat_obj, 
        interactions_db     = liana_db, 
        annot               = params.annot, 
        n_perm              = params.n_perm,
        min_cells           = params.min_cells,
        min_pct             = params.min_pct, 
    )

    POSTPROCESSING_LIANA(
    INFER_LIANA.out.liana_obj, 
        ref_db              = ref_db
    )

    emit: 
    POSTPROCESSING_LIANA.out
}

workflow CELL2CELL {
    take: 
    preprocessing_mtx_dir
    metadata_csv
    cell2cell_db
    ref_db
    
    main: 
    INFER_CELL2CELL(
        preprocessing_mtx_dir,
        meta                = metadata_csv, 
        interactions_db     = cell2cell_db, 
        annot               = params.annot, 
        n_perm              = params.n_perm
    )
    
    POSTPROCESSING_CELL2CELL( 
        INFER_CELL2CELL.out.cell2cell_obj,
        ref_db              = ref_db
    )

    emit: 
    POSTPROCESSING_CELL2CELL.out

}