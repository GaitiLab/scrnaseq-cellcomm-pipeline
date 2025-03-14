include { LIANA_RUN    } from '../../../modules/local/liana/run'
include { LIANA_FORMAT } from '../../../modules/local/liana/format'

workflow LIANA {
    take:
    seurat_obj_prepped
    annot
    n_perm
    min_cells
    min_pct

    main:
    ch_versions = Channel.empty()
    ch_liana_db = Channel.fromPath(params.liana_db)

    ch_ref_db = Channel.fromPath(params.ref_db)
    ch_input = seurat_obj_prepped.combine(ch_liana_db)

    LIANA_RUN(
        ch_input,
        annot,
        n_perm,
        min_cells,
        min_pct,
    )
    ch_versions = ch_versions.mix(LIANA_RUN.out.versions.first())

    LIANA_FORMAT(
        LIANA_RUN.out.rds.combine(ch_ref_db)
    )
    ch_versions = ch_versions.mix(LIANA_FORMAT.out.versions.first())

    emit:
    rds      = LIANA_FORMAT.out.rds
    versions = ch_versions
}
