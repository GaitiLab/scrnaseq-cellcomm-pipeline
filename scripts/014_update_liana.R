# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
cmd_args <- commandArgs(trailingOnly = FALSE)
has_script_filepath <- startsWith(cmd_args, "--file=")
if (sum(has_script_filepath)) {
    setwd(dirname(unlist(strsplit(cmd_args[has_script_filepath], "=")))[2])
}

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Update LIANA database w/ CellChat and CellphoneDB",
    )
    parser$add_argument("--cellchat_db",
        default = "",
        type = "character", help = "Path to CellChat database"
    )
    parser$add_argument("--cpdb_db",
        default = "",
        type = "character", help = "Path to CellphoneDB database"
    )
    parser$add_argument("--liana_db",
        default = "",
        type = "character", help = "Path to LIANA database"
    )

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/001_data_local/interactions_db_v2/")
    args$cellchat_db <- "001_data_local/interactions_db_v2/cellchat_liana_format.rds"
    args$cpdb_db <- "001_data_local/interactions_db_v2/cpdbv5_liana_format.rds"
    args$liana_db <- "001_data_local/interactions_db_v2/liana_db_base.rds"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

args$cellchat_db <- ifelse(file.exists(args$cellchat_db), args$cellchat_db, glue("{args$output_dir}/cellchat_liana_format.rds"))
args$cpdb_db <- ifelse(file.exists(args$cpdb_db), args$cpdb_db, glue("{args$output_dir}/cpdbv5_liana_format.rds"))
args$liana_db <- ifelse(file.exists(args$liana_db), args$liana_db, glue("{args$output_dir}/liana_db_base.rds"))

log_info("Load databases: CellChat and CellphoneDB...")
cellchat_db <- readRDS(args$cellchat_db) %>% mutate(method = "CellChat extracted")
cpdb_db <- readRDS(args$cpdb_db) %>% mutate(method = "CellphoneDB extracted")

log_info("Load LIANA database (Consensus + Ramilowski 2015)...")
liana_db_base <- readRDS(args$liana_db)

log_info("Combine CellChat and CellphoneDB databases...")
cpdb_cellchat_db <- rbind(cpdb_db, cellchat_db)

missing_cols <- setdiff(colnames(liana_db_base), colnames(cpdb_cellchat_db))
log_info(glue("Number of missing columns: {length(missing_cols)}"))
cpdb_cellchat_db[missing_cols] <- ""

log_info("Combine LIANA database with CellChat and CellphoneDB databases...")
liana_db_updated <- rbind(liana_db_base, cpdb_cellchat_db)

log_info(glue("Number of interactions in LIANA database before update: {nrow(liana_db_base)}"))
log_info(glue("Number of interactions in LIANA database after update: {nrow(liana_db_updated)}"))

liana_db_updated <- liana_db_updated %>%
    separate(source_genesymbol, paste0("ligand_subunit_", seq_len(5)), sep = "_", remove = FALSE) %>%
    separate(target_genesymbol, paste0("receptor_subunit_", seq_len(5)), sep = "_", remove = FALSE) %>%
    # define simple interaction (no subunits involved)
    # define complex interaction (subunits involved)
    mutate(
        simple_interaction = paste0(ligand_subunit_1, "_", receptor_subunit_1),
        complex_interaction = paste0(source_genesymbol, "_", target_genesymbol)
    )
# Replace NA with empty string
liana_db_updated[is.na(liana_db_updated)] <- ""

# log_info("Save LIANA database...")
# saveRDS(liana_db_updated, glue("{args$output_dir}/liana_db_combined.rds"))


# ---- Dealing with duplicates ---- #
# If interaction A-B also has the same interaction, but with the specified subunits keep the one with subunits -> more accurate/reflects biology
# Get all interactions that are complex
complex_interactions <- liana_db_updated %>%
    filter(simple_interaction != complex_interaction)
complex_interactions_unique <- complex_interactions %>%
    distinct(complex_interaction, .keep_all = TRUE)
log_info(glue("Number of complex interactions: {nrow(complex_interactions)}"))
log_info(glue("Number of complex interactions after removing duplicates: {nrow(complex_interactions_unique)}"))

# simple_interaction == complex_interaction -> true simple interaction (no subunits involved)
simple_interactions <- liana_db_updated %>%
    filter(simple_interaction == complex_interaction)
simple_interactions_unique <- simple_interactions %>%
    distinct(simple_interaction, .keep_all = TRUE)
log_info(glue("Number of simple interactions: {nrow(simple_interactions)}"))
log_info(glue("Number of simple interactions after removing duplicates: {nrow(simple_interactions_unique)}"))

liana_db_filtered <- rbind(complex_interactions_unique, simple_interactions_unique)
liana_db_filtered <- liana_db_filtered %>% distinct(complex_interaction, .keep_all = TRUE)

log_info(glue("Number of interactions to add: {nrow(liana_db_filtered)}"))
# Should be the same number
log_info(glue("Number of interactions to add after removing duplicates: {nrow(liana_db_filtered)}"))

log_info("Save LIANA database...")
saveRDS(liana_db_filtered, glue("{args$output_dir}/liana_db_combined.rds"))

log_info("COMPLETED!")
