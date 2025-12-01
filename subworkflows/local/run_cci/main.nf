include { CELL2CELL } from '../cell2cell/main.nf'
include { CELLCHAT  } from '../cellchat/main.nf'
include { CPDB      } from '../cpdb/main.nf'
include { LIANA     } from '../liana/main.nf'
workflow RUN_CCI {
    take:
    mtx_dir_prepped
    seurat_obj_prepped
    metadata_csv

    main:
    ch_versions = channel.empty()

    // Set channels for cci 
    ch_cell2cell = channel.empty()
    ch_cellchat = channel.empty()
    ch_cellphonedb = channel.empty()
    ch_liana = channel.empty()

    def cci_tools = params.cci_tools ? params.cci_tools.split(',').collect { it.trim().toLowerCase() } : []

    if (cci_tools.contains("cell2cell")) {
        CELL2CELL(
            mtx_dir_prepped,
            metadata_csv,
            params.annot,
            params.n_perm,
        )
        ch_versions = ch_versions.mix(CELL2CELL.out.versions)
        ch_cell2cell = CELL2CELL.out.rds
    }

    if (cci_tools.contains("cellchat")) {
        CELLCHAT(
            seurat_obj_prepped,
            params.annot,
            params.n_perm,
            params.min_cells,
        )
        ch_versions = ch_versions.mix(CELLCHAT.out.versions)
        ch_cellchat = CELLCHAT.out.rds
    }

    if (cci_tools.contains("cellphonedb")) {
        CPDB(
            mtx_dir_prepped,
            metadata_csv,
            params.annot,
            params.n_perm,
            params.min_pct,
        )
        ch_versions = ch_versions.mix(CPDB.out.versions)
        ch_cellphonedb = CPDB.out.rds
    }
    if (cci_tools.contains("liana")) {
        LIANA(
            seurat_obj_prepped,
            params.annot,
            params.n_perm,
            params.min_cells,
            params.min_pct,
        )

        ch_versions = ch_versions.mix(LIANA.out.versions)
        ch_liana = LIANA.out.rds
    }

    emit:
    cci      = ch_cell2cell.join(ch_cellchat).join(ch_cellphonedb).join(ch_liana)
    versions = ch_versions
}
