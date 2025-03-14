include { CPDB_RUN    } from "../../../modules/local/cpdb/run"
include { CPDB_FORMAT } from "../../../modules/local/cpdb/format"


workflow CPDB {
    take:
    mtx_dir_prepped
    metadata_csv
    annot
    n_perm
    min_pct

    main:
    ch_versions = Channel.empty()
    ch_cellphone_db = Channel.fromPath(file(params.cellphonedb_db))
    ch_ref_db = Channel.fromPath(file(params.ref_db))

    ch_input = mtx_dir_prepped.combine(metadata_csv).combine(ch_cellphone_db)
    CPDB_RUN(
        ch_input,
        annot,
        n_perm,
        min_pct,
    )
    ch_versions = ch_versions.mix(CPDB_RUN.out.versions.first())

    CPDB_FORMAT(
        CPDB_RUN.out.txt.combine(ch_ref_db)
    )
    ch_versions = ch_versions.mix(CPDB_FORMAT.out.versions.first())

    emit:
    rds      = CPDB_FORMAT.out.rds
    versions = ch_versions
}
