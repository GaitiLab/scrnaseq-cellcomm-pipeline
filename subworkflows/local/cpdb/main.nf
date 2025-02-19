include { RUN_CPDB    } from "../../../modules/local/runcpdb.nf"
include { FORMAT_CPDB } from "../../../modules/local/formatcpdb.nf"


workflow CPDB {
    take:
    mtx_dir_prepped
    metadata_csv   
    cellphone_db   
    ref_db         
    annot          
    n_perm         
    min_pct        

    main:
    RUN_CPDB(
        mtx_dir_prepped,
        metadata_csv,
        cellphone_db,
        annot,
        n_perm,
        min_pct
    )

    FORMAT_CPDB(
        RUN_CPDB.out.txt,
        ref_db
    )

    emit:
    rds = FORMAT_CPDB.out.rds
}
