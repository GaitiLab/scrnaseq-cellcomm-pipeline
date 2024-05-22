include { formatting_cellchat; formatting_liana; formatting_cell2cell; formatting_cpdb } from "../nf-modules/formatting_cci.nf"

workflow POSTPROCESSING_CCI {
    take: 
    liana_obj
    cell2cell_obj
    cpdb_obj
    cellchat_obj
    ref_db 

    main: 
    formatting_cellchat(
        cellchat_obj, 
        interactions_db     = ref_db
    )

    formatting_liana(
        liana_obj, 
        ref_db              = ref_db
    )

    formatting_cell2cell(
        cell2cell_obj, 
        ref_db              = ref_db
    )

    formatting_cpdb(
        cpdb_obj, 
        interactions_db     = ref_db
    )

    emit: 
    cci_liana       = formatting_liana.out 
    cci_cell2cell   = formatting_cell2cell.out
    cci_cpdb        = formatting_cpdb.out 
    cci_cellchat    = formatting_cellchat.out
}