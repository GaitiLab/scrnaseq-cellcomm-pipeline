include { CELLCHAT_RUN    } from "../../../modules/local/cellchat/run"
include { CELLCHAT_FORMAT } from "../../../modules/local/cellchat/format"

workflow CELLCHAT {
    take:
    seurat_obj_prepped
    annot
    n_perm
    min_cells

    main:
    ch_versions = Channel.empty()
    ch_cellchat_db = Channel.fromPath(params.cellchat_db)
    ch_ref_db = Channel.fromPath(params.ref_db)

    ch_input = seurat_obj_prepped.combine(ch_cellchat_db)

    CELLCHAT_RUN(
        ch_input,
        annot,
        n_perm,
        min_cells,
    )
    ch_versions = ch_versions.mix(CELLCHAT_RUN.out.versions.first())

    CELLCHAT_FORMAT(
        CELLCHAT_RUN.out.rds.combine(ch_ref_db)
    )
    ch_versions = ch_versions.mix(CELLCHAT_FORMAT.out.versions.first())

    emit:
    rds      = CELLCHAT_FORMAT.out.rds
    versions = ch_versions
}
