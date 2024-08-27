include { infer_cellchat; infer_liana; infer_cell2cell; infer_cpdb } from "../nf-modules/infer_interactions.nf"
include { formatting_cellchat; formatting_liana; formatting_cell2cell; formatting_cpdb } from "../nf-modules/formatting_cci.nf"

workflow CELLCHAT {
    take: 
    preprocessing_seurat_obj
    cellchat_db
    ref_db 

    main: 
    infer_cellchat(
        preprocessing_seurat_obj, 
        interactions_db     = cellchat_db, 
        annot               = params.annot,
        n_perm              = params.n_perm,
        min_cells           = params.min_cells
    )

    formatting_cellchat(
        infer_cellchat.out.cellchat_obj, 
        interactions_db     = ref_db
    )

    emit: 
    formatting_cellchat.out
}

workflow CPDB {
    take: 
    preprocessing_mtx_dir
    metadata_csv
    cellphone_db
    ref_db
    
    main: 
    infer_cpdb(
        preprocessing_mtx_dir,
        meta                = metadata_csv, 
        interactions_db     = cellphone_db, 
        annot               = params.annot, 
        n_perm              = params.n_perm, 
        min_pct             = params.min_pct, 
    )

    formatting_cpdb(
        infer_cpdb.out.cpdb_obj, 
        interactions_db     = ref_db
    )

    emit: 
    formatting_cpdb.out
}

workflow LIANA {
    take: 
    preprocessing_seurat_obj
    liana_db
    ref_db
    
    main: 
    infer_liana(
        preprocessing_seurat_obj, 
        interactions_db     = liana_db, 
        annot               = params.annot, 
        n_perm              = params.n_perm,
        min_cells           = params.min_cells,
        min_pct             = params.min_pct, 
    )

    formatting_liana(
        infer_liana.out.liana_obj, 
        ref_db              = ref_db
    )

    emit: 
    formatting_liana.out
}

workflow CELL2CELL {
    take: 
    preprocessing_mtx_dir
    metadata_csv
    cell2cell_db
    ref_db
    
    main: 
    infer_cell2cell(
        preprocessing_mtx_dir,
        meta                = metadata_csv, 
        interactions_db     = cell2cell_db, 
        annot               = params.annot, 
        n_perm              = params.n_perm
    )
    
    formatting_cell2cell( 
        infer_cell2cell.out.cell2cell_obj,
        ref_db              = ref_db
    )

    emit: 
    formatting_cell2cell.out

}