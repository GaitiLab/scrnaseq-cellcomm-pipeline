include { RUN_LIANA } from "../../../modules/local/runliana.nf"
include { FORMAT_LIANA } from "../../../modules/local/formatliana.nf"

workflow LIANA {
    take:
    seurat_obj_prepped
    liana_db
    ref_db
    annot
    n_perm
    min_cells
    min_pct

    main:
    RUN_LIANA(
        seurat_obj_prepped,
        liana_db,
        annot,
        n_perm,
        min_cells,
        min_pct
    )

    FORMAT_LIANA(
        RUN_LIANA.out.rds,
        ref_db
    )

    emit:
    rds = FORMAT_LIANA.out.rds
}
