include { POSTPROCESSING_CELLCHAT; POSTPROCESSING_LIANA; POSTPROCESSING_CELL2CELL; POSTPROCESSING_CPDB } from "../nf-modules/filtering.nf"

workflow POSTPROCESSING_CCI {
    take: 
    liana_obj
    cell2cell_obj
    cpdb_obj
    cellchat_obj
    ref_db 

    main: 
    POSTPROCESSING_CELLCHAT(
        cellchat_obj, 
        interactions_db     = ref_db
    )

    POSTPROCESSING_LIANA(
        liana_obj, 
        ref_db              = ref_db
    )

    POSTPROCESSING_CELL2CELL(
        cell2cell_obj, 
        ref_db              = ref_db
    )

    POSTPROCESSING_CPDB(
        cpdb_obj, 
        interactions_db     = ref_db
    )

    emit: 
    cci_liana       = POSTPROCESSING_LIANA.out 
    cci_cell2cell   = POSTPROCESSING_CELL2CELL.out
    cci_cpdb        = POSTPROCESSING_CPDB.out 
    cci_cellchat    = POSTPROCESSING_CELLCHAT.out
}