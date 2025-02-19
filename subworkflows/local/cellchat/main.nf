include { RUN_CELLCHAT } from "../../../modules/local/runcellchat.nf"
include { FORMAT_CELLCHAT } from "../../../modules/local/formatcellchat.nf"

workflow CELLCHAT {
    take:
    seurat_obj_prepped
    cellchat_db
    ref_db
    annot
    n_perm
    min_cells

    main:
    RUN_CELLCHAT(
        seurat_obj_prepped,
        cellchat_db,
        annot,
        n_perm,
        min_cells
    )

    FORMAT_CELLCHAT(
        RUN_CELLCHAT.out.rds,
        ref_db
    )

    emit:
    rds = FORMAT_CELLCHAT.out.rds
}
