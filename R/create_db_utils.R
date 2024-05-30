#' @title Check genesymbol length
#' @param vec vector with genesymbols
#' @export
check_length_genesymbol <- function(vec) {
    mask <- (nchar(vec[!is.na(vec)]) > 1)
    sum(mask) == length(mask)
}

#' @title order vector
#' @param vec vector
#' @description order elements in vector alphabetically and return vector with same length
#' @export
order_vector <- function(vec) {
    sorted_vec <- sort(vec)
    placeholder <- rep(NA, length(vec))
    placeholder[seq(length(sorted_vec))] <- sorted_vec
    return(placeholder)
}

#' @title make complex
#' @param vec vector
#' @description Merge items in vector by collapse character specified by user, default = ':'
#' @export
make_complex <- function(vec, collapse = ":") {
    paste0(vec[!is.na(vec)], collapse = collapse)
}

#' @title Look up uniprot based on genesymbol
#' @param gene genesymbol
#' @param gene_uniprot_table dataframe containing at least 'genesymbol' and 'uniprot' columns
#' @export
#' @importFrom dplyr %>%
lookup_uniprot <- function(gene, gene_uniprot_table) {
    gene_uniprot_table %>%
        filter(genesymbol == gene) %>%
        dplyr::select(uniprot)
}

#' @title Check mapping
#' @description Check whether all entries in vector are found in a dataframe
#' @param vec vector
#' @param lookup_table dataframe for looking up valuees
#' @param lookup_var variable in lookup_table for checking overlap with vec
#' @return returns whether all entries of the vector where found in the lookup_table (boolean)
is_mapped <- function(vec, lookup_table, lookup_var = "genesymbol") {
    vec_not_na <- vec[!is.na(vec)]
    # Number of entries found should match with the length of vector without missing values
    return(
        length(
            # Check number of entries are found in table
            intersect(
                vec_not_na,
                lookup_table %>% dplyr::pull(!!dplyr::sym(lookup_var))
            )
        ) == length(vec_not_na)
    )
}

#' @title Add proteins
#' @param vec vector
#' @param gene_uniprot_table dataframe containing at least 'genesymbol' and 'uniprot' columns
#' @return string with all protein complexes that are involved separate by an underscore, if no protein found then returns vector of NAs (same length as input vector)
add_proteins <- function(vec, gene_uniprot_table) {
    tmp <- paste0(lapply(vec[!is.na(vec)], lookup_uniprot, gene_uniprot_table = gene_uniprot_table) %>% unlist(), collapse = "_")
    if (tmp == "") {
        return(NA)
    } else {
        tmp
    }
}

#' @title Load LIANA DB
#' @description Loads and combines the LIANA Consensus and Ramilowski 2015 database
#' @return dataframe with the two databases
#' @export
#' @importFrom dplyr %>%
load_liana_dbs <- function() {
    # Loading Consensus database
    liana_consensus <- liana::select_resource("Consensus")[[1]] %>% dplyr::mutate(method = "LIANA consensus")

    # load Ramilowski 2015 database
    liana_ramilowski <- liana::select_resource("Ramilowski2015")[[1]] %>% dplyr::mutate(method = "LIANA Ramilowski2015")

    # Account for missing columns
    missing_liana <- setdiff(colnames(liana_consensus), colnames(liana_ramilowski))
    missing_ramilowski <- setdiff(colnames(liana_ramilowski), colnames(liana_consensus))

    liana_consensus[missing_ramilowski] <- ""
    liana_ramilowski[missing_liana] <- ""

    # Combining the two LIANA databases
    liana_db <- rbind(liana_consensus, liana_ramilowski)

    return(liana_db)
}

#' @title Load CellChat DB
#' @description load cellcat database and format for merging
#' @return dataframe with 'source_genesymbol' and 'target_genesymbol'
#' @importFrom dplyr %>%
#' @export
load_cellchat_db <- function() {
    cellchat_db <- CellChat::CellChatDB.human
    cellchat_db_interaction <- cellchat_db$interaction

    interactions <- cellchat_db_interaction %>%
        tibble::remove_rownames() %>%
        dplyr::select(interaction_name, ligand, receptor, ligand.symbol, receptor.symbol) %>%
        dplyr::mutate(
            # Separate subunits with underscores '_'
            source_genesymbol = stringr::str_replace_all(ligand.symbol, ", ", "_"),
            target_genesymbol = stringr::str_replace_all(receptor.symbol, ", ", "_")
        ) %>%
        dplyr::mutate(
            # Ensure all genes are capitalized
            source_genesymbol = ifelse(source_genesymbol == "", toupper(ligand), source_genesymbol),
            target_genesymbol = ifelse(target_genesymbol == "", toupper(receptor), target_genesymbol)
        ) %>%
        dplyr::select(source_genesymbol, target_genesymbol)
    return(interactions)
}

#' @title Load CellPhoneDB DB
#' @description load CellPhoneDB database and format for merging
#' @param source_cpdb_dir path to 'interaction_input.csv' from CellPhoneDB database
#' @return dataframe with 'source_genesymbol' and 'target_genesymbol'
#' @importFrom dplyr %>%
#' @export
load_cpdb_db <- function(source_cpdb_dir) {
    cpdb_interaction_input <- read.csv(
        glue::glue("{source_cpdb_dir}/interaction_input.csv"),
        sep = ","
    )
    # Split subunits by '-'
    interaction_units <- stringr::str_split(cpdb_interaction_input$interactors, "-", simplify = TRUE)
    # Add column names for subunits
    colnames(interaction_units) <- paste0("unit_", seq_len(ncol(interaction_units)))

    custom_cellphonedb <- data.frame(interaction_units) %>%
        dplyr::mutate(
            # Handle individual genes with a '-', i.e. merging with previous subunit if length of unit_2 is a single character
            source_genesymbol_complex = dplyr::case_when(
                nchar(unit_2) == 1 ~ paste0(unit_1, "-", unit_2), TRUE ~ unit_1
            ),
            target_genesymbol_complex = dplyr::case_when(
                nchar(unit_2) == 1 ~ unit_3,
                nchar(unit_3) == 1 ~ paste0(unit_2, "-", unit_3), TRUE ~ unit_2
            )
        ) %>%
        # Format source_genesymbol and target_genesymbol
        dplyr::mutate(
            source_genesymbol = stringr::str_replace_all(source_genesymbol_complex, "\\+", "_"),
            target_genesymbol = stringr::str_replace_all(target_genesymbol_complex, "\\+", "_")
        ) %>%
        dplyr::select(source_genesymbol, target_genesymbol)
    return(custom_cellphonedb)
}

#' @title Hard combine databases
#' @description Combine databases based on 'source_genesymbol' and 'target_genesymbol'
#' @param cellchat_db cellchat database
#' @param cpdb_db cpdb_database
#' @param liana_db liana database
#' @export
#' @importFrom dplyr %>%
hard_combine_dbs <- function(cellchat_db, cpdb_db, liana_db) {
    # Load databases: CellChat and CellphoneDB
    cellchat_db <- cellchat_db %>% dplyr::mutate(method = "CellChat extracted")
    cpdb_db <- cpdb_db %>% dplyr::mutate(method = "CellphoneDB extracted")

    # Combine CellChat and CellphoneDB databases
    cpdb_cellchat_db <- rbind(cpdb_db, cellchat_db)

    # Add missing columns and set to empty string
    missing_cols <- setdiff(colnames(liana_db), colnames(cpdb_cellchat_db))
    cpdb_cellchat_db[missing_cols] <- ""

    # Combine LIANA database with CellChat and CellphoneDB databases
    db_combined <- rbind(liana_db, cpdb_cellchat_db) %>%
        dplyr::mutate(
            complex_interaction = paste0(source_genesymbol, "_", target_genesymbol)
        )
    return(db_combined)
}

#' @title Filter and clean-up database
#' @param db dataframe with database of interactions, the following columns should at least be present 'source_genesymbol' and 'target_genesymbol'
#' @export
#' @importFrom dplyr %>%
db_filtering <- function(db) {
    db <- db %>%
        dplyr::distinct(complex_interaction, .keep_all = TRUE) %>%
        # Split into subunits
        tidyr::separate(source_genesymbol, paste0("source_genesymbol_subunit_", seq_len(5)),
            sep = "_", remove = FALSE
        ) %>%
        tidyr::separate(target_genesymbol, paste0("target_genesymbol_subunit_", seq_len(5)),
            sep = "_", remove = FALSE
        )
    # Set empty strings to NA
    db[db == ""] <- NA

    # Order subunits alphabetically and replace
    source_genesymbol_subunits_ordered <- data.frame(t(apply(db %>% dplyr::select(dplyr::starts_with("source_genesymbol_subunit_")), 1, order_vector)))
    target_genesymbol_subunits_ordered <- data.frame(t(apply(db %>% dplyr::select(dplyr::starts_with("target_genesymbol_subunit_")), 1, order_vector)))

    colnames(source_genesymbol_subunits_ordered) <- paste0("source_genesymbol_subunit_", seq(5))
    colnames(target_genesymbol_subunits_ordered) <- paste0("target_genesymbol_subunit_", seq(5))

    # Make complexes based on ordered subunits
    source_genesymbol_complex <- apply(source_genesymbol_subunits_ordered, 1, make_complex)
    target_genesymbol_complex <- apply(target_genesymbol_subunits_ordered, 1, make_complex)

    # Add ordered subunits to dataframe (replaces)
    db <- cbind(
        db %>% dplyr::select(
            -dplyr::all_of(c(
                paste0("source_genesymbol_subunit_", seq(5)),
                paste0("target_genesymbol_subunit_", seq(5))
            ))
        ),
        cbind(source_genesymbol_subunits_ordered, target_genesymbol_subunits_ordered)
    )

    db <- db %>%
        # Keep unordered (original) 'source_genesymbol' and 'target_genesymbol'
        dplyr::rename(
            source_genesymbol_old = source_genesymbol,
            target_genesymbol_old = target_genesymbol
        ) %>%
        dplyr::mutate(
            source_genesymbol = source_genesymbol_complex,
            target_genesymbol = target_genesymbol_complex
        ) %>%
        # Remove duplicates based on ordered complexes
        dplyr::distinct(source_genesymbol, target_genesymbol, .keep_all = TRUE)

    # Check if subunits have at least 2 characters (otherwise invalid gene)
    target_genesymbol_mask <- apply(db %>% dplyr::select(dplyr::starts_with("target_genesymbol_subunit_")), 1, check_length_genesymbol)
    source_genesymbol_mask <- apply(db %>% dplyr::select(dplyr::starts_with("source_genesymbol_subunit_")), 1, check_length_genesymbol)
    # Remove interactions where genes only have a single character
    db <- db[source_genesymbol_mask & target_genesymbol_mask, ]

    # Adding protein ids (uniprot)
    all_genes <- db %>%
        dplyr::select(source_genesymbol, target_genesymbol) %>%
        as.list() %>%
        unlist() %>%
        stringr::str_split(., ":") %>%
        unlist() %>%
        unique()
    all_genes <- data.frame(genesymbol = all_genes[all_genes != ""] %>% unique())

    # CellPhoneDB needs: gene_name, uniprot, hgnc_symbol, ensemble, uniprot, protein_name
    proteins <- OmnipathR::translate_ids(
        # Query
        all_genes,
        # Match ID
        genesymbol = genesymbol,
        # Variables to retrieve from OmnipathR
        uniprot, uniprot_entry
    )
    proteins <- OmnipathR::translate_ids(
        # Query
        proteins,
        # Match ID
        uniprot = uniprot,
        ensembl
    )
    # Remove duplicates based on genesymbol
    proteins <- proteins %>%
        na.omit() %>%
        dplyr::distinct(genesymbol, .keep_all = TRUE)

    # Check mapping gene to protein (uniprot)
    source_mask <- apply(
        db %>% dplyr::select(dplyr::starts_with("source_genesymbol_subunit_")), 1, is_mapped,
        lookup_table = proteins,
        lookup_var = "genesymbol"
    )

    target_mask <- apply(db %>% dplyr::select(dplyr::starts_with("target_genesymbol_subunit_")), 1,
        is_mapped,
        lookup_table = proteins,
        lookup_var = "genesymbol"
    )
    # Remove interactions for which we don't have a matching uniprot (protein)
    db <- db[source_mask & target_mask, ]

    # Add proteins to dataframe
    db$target <- apply(
        db %>% dplyr::select(dplyr::starts_with("target_genesymbol_subunit_")),
        1,
        add_proteins,
        gene_uniprot_table = proteins
    )

    db$source <- apply(
        db %>% dplyr::select(dplyr::starts_with("source_genesymbol_subunit_")),
        1,
        add_proteins,
        gene_uniprot_table = proteins
    )

    # Ensure that the dataframe only contains interactions for which the protein IDs are available for both source AND target
    db <- db %>%
        filter(!is.na(source), !is.na(target)) %>%
        # Rename previously set 'complex_interaction' (unordered subunits)
        rename(complex_interaction_OLD = complex_interaction) %>%
        dplyr::mutate(
            # Separate complexes by underscores '_'
            source_genesymbol = stringr::str_replace_all(source_genesymbol, "\\:", "_"),
            target_genesymbol = stringr::str_replace_all(target_genesymbol, "\\:", "_"),
            complex_interaction = paste0(source_genesymbol, "__", target_genesymbol)
        )
    return(list(db = db, gene_info_table = proteins))
}

#' @title Update CellPhoneDB
#' @description Update CellPhoneDB and generate required files
#' @param source_cpdb_dir directory with CellPhoneDB files
#' @param db database with at least 'source_genesymbol' and 'target_genesymbol' and their subunits starting with 'source_genesymbol_subunit' or 'target_genesymbol_subunit'
#' @param gene_info_table dataframe with uniprot, protein name, genesymbol and ensembl
#' @param output_dir path to directory for saving the generated CellPhoneDB files
#' @param return_list return list with generated dataframes (default = FALSE)
#' @return if return_list = TRUE, return list with complex_input, gene_input, interaction_input and protein_input
#' @export
#' @importFrom dplyr %>%
update_cpdb_db <- function(source_cpdb_dir, db, gene_info_table, output_dir, return_list = FALSE) {
    # Grab column names for generating new files
    complex_input_cols <- colnames(read.csv(glue::glue("{source_cpdb_dir}/complex_input.csv")))
    gene_input_cols <- colnames(read.csv(glue::glue("{source_cpdb_dir}/gene_input.csv")))
    interaction_input_cols <- colnames(read.csv(glue::glue("{source_cpdb_dir}/interaction_input.csv")))
    protein_input_cols <- colnames(read.csv(glue::glue("{source_cpdb_dir}/protein_input.csv")))

    gene_info_table <- gene_info_table %>% dplyr::distinct(uniprot, .keep_all = TRUE)

    # Gene input
    # mandatory: gene_name, uniprot, hgnc_symbol and ensemble
    genes_input <- gene_info_table %>%
        dplyr::rename(hgnc_symbol = genesymbol, ensemble = ensembl) %>%
        dplyr::mutate(
            gene_name = hgnc_symbol
        ) %>%
        dplyr::select(
            gene_name, uniprot, hgnc_symbol, ensemble
        )

    # Protein input
    # mandatory: uniprot, protein_name
    protein_input <- gene_info_table %>%
        dplyr::rename(protein_name = "uniprot_entry") %>%
        dplyr::select(uniprot, protein_name) %>%
        dplyr::distinct(uniprot, .keep_all = TRUE)

    # Complex input
    # mandatory: complex_name, complex_subunits
    ligand_complexes <- db %>%
        tidyr::separate(source, paste0("uniprot_", seq_len(5)), sep = "_", remove = FALSE) %>%
        dplyr::select(source_genesymbol, dplyr::starts_with("uniprot_")) %>%
        dplyr::rename(complex_name = source_genesymbol)
    receptor_complexes <- db %>%
        tidyr::separate(target, paste0("uniprot_", seq_len(5)), sep = "_", remove = FALSE) %>%
        dplyr::select(target_genesymbol, dplyr::starts_with("uniprot_")) %>%
        dplyr::rename(complex_name = target_genesymbol)
    complex_input <- rbind(ligand_complexes, receptor_complexes) %>% dplyr::distinct()
    complex_input[complex_input == ""] <- NA

    # Ensure that the dataframe only contains protein complexes, i.e. uniprot_2 should NOT be empty
    complex_input <- complex_input %>% filter(!is.na(uniprot_2))

    # Order subunits from complexes
    complex_subunits_ordered <- data.frame(t(apply(complex_input %>% dplyr::select(dplyr::starts_with("uniprot_")), 1, order_vector)))
    colnames(complex_subunits_ordered) <- paste0("uniprot_", seq_len(5))

    # Create complexes, merge by '_'
    dummy <- apply(complex_subunits_ordered, 1, make_complex, collapse = "_")

    # Only keep unique complexes (based on ordered subunits)
    complex_input <- cbind(complex_input %>% dplyr::select(-dplyr::all_of(paste0("uniprot_", seq_len(5)))), complex_subunits_ordered) %>%
        dplyr::mutate(dummy = dummy) %>%
        dplyr::distinct(dummy, .keep_all = TRUE) %>%
        dplyr::select(-dummy)

    # Interaction Input
    # mandatory:  “partner_a”; “partner_b”; “annotation_strategy”; “source”
    interactions_input <- db %>%
        dplyr::select(
            # Proteins
            source, target,
            # Necessary for checks
            source_genesymbol_subunit_2, target_genesymbol_subunit_2,
            # Genesymbols
            source_genesymbol, target_genesymbol
        ) %>%
        # Rename protein columns
        dplyr::rename(partner_a = source, partner_b = target) %>%
        dplyr::mutate(
            annotation_strategy = "user_curated",
            version = "CellPhoneDBcore4.1",
            sources = "User curated",
            # If complex, then user 'complex_name', i.e. source_genesymbol/target_genesymbol
            partner_a = ifelse(is.na(source_genesymbol_subunit_2), partner_a, source_genesymbol),
            partner_b = ifelse(is.na(target_genesymbol_subunit_2), partner_b, target_genesymbol),
        ) %>%
        # Select
        dplyr::select(
            partner_a, partner_b, annotation_strategy, sources
        )

    # account for missing cols
    missing_cols <- setdiff(gene_input_cols, colnames(genes_input))
    genes_input[missing_cols] <- ""
    missing_cols <- setdiff(protein_input_cols, colnames(protein_input))
    protein_input[missing_cols] <- ""
    missing_cols <- setdiff(complex_input_cols, colnames(complex_input))
    complex_input[missing_cols] <- ""
    complex_input$version <- "CellPhoneDBcore4.1"
    missing_cols <- setdiff(interaction_input_cols, colnames(interactions_input))
    interactions_input[missing_cols] <- ""

    genes_input[is.na(genes_input)] <- ""
    protein_input[is.na(protein_input)] <- ""
    complex_input[is.na(complex_input)] <- ""
    interactions_input[is.na(interactions_input)] <- ""

    # Save files
    write.csv(genes_input, glue::glue("{output_dir}/gene_input.csv"), row.names = FALSE, quote = FALSE)
    write.csv(protein_input, glue::glue("{output_dir}/protein_input.csv"), row.names = FALSE, quote = FALSE)
    write.csv(complex_input, glue::glue("{output_dir}/complex_input.csv"), row.names = FALSE, quote = FALSE)
    write.csv(interactions_input, glue::glue("{output_dir}/interaction_input.csv"), row.names = FALSE, quote = FALSE)

    if (return_list) {
        return(list(genes_input = genes_input, protein_input = protein_input, complex_input = complex_input, interactions_input = interactions_input))
    }
}

#' @title Update CellChat DB using (updated) CellPhoneDB
#' @param cpdb_dir directory with CellPhoneDB files
#' @param output_dir path to directory for saving the generated CellPhoneDB database (.rds)
#' @param db database
#' @export
#' @importFrom dplyr %>%
update_cellchat_db <- function(cpdb_dir, db, output_dir) {
    # Ref: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Update-CellChatDB.html
    cellchat_db <- CellChat::CellChatDB.human
    cellchat_db_interaction <- cellchat_db$interaction
    cellchat_db_complex <- cellchat_db$complex
    cellchat_db_cofactor <- cellchat_db$cofactor
    cellchat_db_geneinfo <- cellchat_db$geneInfo

    # Load CellPhoneDB files
    geneInfo <- read.csv(file = glue::glue("{cpdb_dir}/gene_input.csv"))
    geneInfo$Symbol <- geneInfo$hgnc_symbol
    geneInfo <- dplyr::select(geneInfo, -c("ensembl"))
    geneInfo <- unique(geneInfo)

    geneInfo <- geneInfo %>%
        dplyr::rename(EntryID.uniprot = uniprot) %>%
        dplyr::select(EntryID.uniprot, Symbol) %>%
        dplyr::mutate()

    # Get missing columns
    missing_cols <- setdiff(colnames(cellchat_db_geneinfo), colnames(geneInfo))
    geneInfo[, missing_cols] <- ""

    cellchat_db_geneinfo <- rbind(cellchat_db_geneinfo, geneInfo) %>% dplyr::distinct(Symbol, .keep_all = TRUE)

    # Create new CellChatDB
    cellchatDB_Omni <- liana:::cellchat_formatDB(
        ccDB = list(
            complex = cellchat_db_complex,
            geneInfo = cellchat_db_geneinfo,
            interaction = cellchat_db_interaction,
            cofactor = cellchat_db_cofactor
        ),
        op_resource = db,
        exclude_anns = c()
    )

    # Save CellChatDB
    saveRDS(cellchatDB_Omni,
        file = glue::glue("{output_dir}/cellchat_db.rds")
    )
}
