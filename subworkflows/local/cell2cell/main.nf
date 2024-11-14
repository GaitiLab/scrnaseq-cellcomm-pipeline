include { RUN_CELL2CELL } from "../../../modules/local/runcell2cell.nf"
include { FORMAT_CELL2CELL } from "../../../modules/local/formatcell2cell.nf"

workflow CELL2CELL {
    take:
    mtx_dir_prepped
    metadata_csv
    cell2cell_db
    ref_db
    annot
    n_perm

    main:
    RUN_CELL2CELL(
        mtx_dir_prepped,
        metadata_csv,
        cell2cell_db,
        annot,
        n_perm
    )

    FORMAT_CELL2CELL(
        RUN_CELL2CELL.out.csv,
        ref_db
    )

    emit:
    rds = FORMAT_CELL2CELL.out.rds
}
