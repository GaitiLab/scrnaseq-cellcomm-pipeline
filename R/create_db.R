#' @title Create interactions databasee
#' @description Generates database for LIANA, CellChat, Cell2Cell and CellPhoneDB
#' @param source_cpdb_dir path to 'interaction_input.csv' from CellPhoneDB database
#' @param output_dir path to directory for saving the generated CellPhoneDB files
#' @param work_dir path to directory with /jobs and /Python
#' @export
create_db <- function(source_cpdb_dir, output_dir, work_dir = here::here()) {
    stopifnot(file.exists(source_cpdb_dir))
    stopifnot(file.exists(output_dir))
    stopifnot(file.exists(work_dir))
    # source_cpdb_dir <- glue("{here::here()}/000_misc_local/references/cellphonedb_v5.0.0/")
    # output_dir <- "output/wip_interactions_db"
    # Get LIANA (Consensus + Ramilowski)
    liana_db <- load_liana_dbs()

    # Get CellChat interactions
    cellchat_db <- load_cellchat_db()

    # Get CellPhoneDB interactions
    cpdb_db <- load_cpdb_db(source_cpdb_dir)

    # Combine databases
    db <- hard_combine_dbs(liana_db = liana_db, cellchat_db = cellchat_db, cpdb_db = cpdb_db)

    # Filter out duplicates + formatting
    db_filtered <- db_filtering(db)
    db_filtered_df <- db_filtered$db

    gene_info_table <- db_filtered$gene_info_table
    # Update CellPhoneDB database + save
    update_cpdb_db(
        source_cpdb_dir = source_cpdb_dir,
        db = db_filtered_df,
        gene_info_table = gene_info_table,
        output_dir = output_dir,
        return_list = FALSE
    )

    # Update CellChat datbase + save
    update_cellchat_db(
        cpdb_dir = output_dir,
        output_dir = output_dir,
        db = db_filtered_df
    )

    # Save LIANA db and Cell2Cell db
    saveRDS(db_filtered_df, glue::glue("{output_dir}/liana_db.rds"))
    write.csv(db_filtered_df, glue::glue("{output_dir}/cell2cell_db.csv"), row.names = FALSE)

    system(glue::glue("{work_dir}/jobs/000_create_db.sh {work_dir} {output_dir}"))
    cpdb_file <- list.files(output_dir, pattern = ".zip", full.names = TRUE)[1]
    file.rename(cpdb_file, glue::glue("{output_dir}/cpdb.zip"))
}
