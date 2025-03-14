include { CELL2CELL_RUN    } from "../../../modules/local/cell2cell/run"
include { CELL2CELL_FORMAT } from "../../../modules/local/cell2cell/format"

workflow CELL2CELL {
    take:
    mtx_dir_prepped
    metadata_csv
    annot
    n_perm

    main:
    ch_versions = Channel.empty()
    ch_ref_db = Channel.fromPath(file(params.ref_db))
    cell2cell_db = Channel.fromPath(file(params.cell2cell_db))

    ch_input = mtx_dir_prepped.combine(metadata_csv).combine(cell2cell_db)

    CELL2CELL_RUN(
        ch_input,
        annot,
        n_perm,
    )
    ch_versions = ch_versions.mix(CELL2CELL_RUN.out.versions.first())

    CELL2CELL_FORMAT(
        CELL2CELL_RUN.out.csv.combine(ch_ref_db)
    )
    ch_versions = ch_versions.mix(CELL2CELL_FORMAT.out.versions.first())

    emit:
    rds      = CELL2CELL_FORMAT.out.rds
    versions = ch_versions
}
