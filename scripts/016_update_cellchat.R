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
        description = "Update CellChat database w/ LIANA",
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
    args$liana_db <- "001_data_local/interactions_db_v2/liana_db.rds"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Load additional libraries
pacman::p_load_gh("jinworks/CellChat")

# Ref: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Update-CellChatDB.html

log_info("Load CellChat database...")
cellchat_db <- CellChatDB.human
cellchat_db_interaction <- cellchat_db$interaction
cellchat_db_complex <- cellchat_db$complex
cellchat_db_cofactor <- cellchat_db$cofactor
cellchat_db_geneinfo <- cellchat_db$geneInfo

log_info("Load LIANA database...")
liana_db <- readRDS(ifelse(file.exists(args$liana_db),
    args$liana_db, glue("{args$output_dir}/liana_db.rds")
))
log_info(glue("Number of interactions: {nrow(liana_db)}"))

# Create new CellChatDB
log_info("Create new CellChatDB...")
cellchatDB_Omni <- liana:::cellchat_formatDB(
    ccDB = CellChat::CellChatDB.human,
    op_resource = liana_db,
    exclude_anns = c()
)

# new_cellchat_db_interaction <- cellchatDB_Omni$interaction
# new_cellchat_db_complex <- cellchatDB_Omni$complex
# new_cellchat_db_cofactor <- cellchatDB_Omni$cofactor
# new_cellchat_db_geneinfo <- cellchatDB_Omni$geneInfo

# head(new_cellchat_db_interaction)


# # Account for missing columns
# missing_cols <- setdiff(colnames(cellchat_db_interaction), colnames(new_cellchat_db_interaction))
# new_cellchat_db_interaction[missing_cols] <- ""
# missing_cols_old <- setdiff(colnames(new_cellchat_db_interaction), colnames(cellchat_db_interaction))
# cellchat_db_interaction[missing_cols_old] <- ""

# missing_cols <- setdiff(colnames(cellchat_db_complex), colnames(new_cellchat_db_complex))
# new_cellchat_db_complex[missing_cols] <- ""
# missing_cols_old <- setdiff(colnames(new_cellchat_db_complex), colnames(cellchat_db_complex))
# cellchat_db_complex[missing_cols_old] <- ""

# missing_cols <- setdiff(colnames(cellchat_db_cofactor), colnames(new_cellchat_db_cofactor))
# new_cellchat_db_cofactor[missing_cols] <- ""
# missing_cols_old <- setdiff(colnames(new_cellchat_db_cofactor), colnames(cellchat_db_cofactor))
# cellchat_db_cofactor[missing_cols_old] <- ""

# missing_cols <- setdiff(colnames(cellchat_db_geneinfo), colnames(new_cellchat_db_geneinfo))
# new_cellchat_db_geneinfo[missing_cols] <- ""
# missing_cols_old <- setdiff(colnames(new_cellchat_db_geneinfo), colnames(cellchat_db_geneinfo))
# cellchat_db_geneinfo[missing_cols_old] <- ""

# # Add new interactions
# cellchat_db_interaction_updated <- rbind(cellchat_db_interaction, new_cellchat_db_interaction)
# log_info(glue("Number of interactions: {nrow(cellchat_db_interaction)}"))
# log_info(glue("Number of interactions combined: {nrow(cellchat_db_interaction_updated)}"))
# cellchat_db_interaction_updated <- cellchat_db_interaction_updated %>%
#     mutate(
#         source_genesymbol = str_replace_all(ligand.symbol, ", ", "_"),
#         target_genesymbol = str_replace_all(receptor.symbol, ", ", "_")
#     ) %>%
#     mutate(
#         source_genesymbol = ifelse(source_genesymbol == "", ligand, source_genesymbol),
#         target_genesymbol = ifelse(target_genesymbol == "", receptor, target_genesymbol)
#     ) %>%
#     mutate(complex_interaction = paste0(source_genesymbol, "_", target_genesymbol))
# cellchat_db_interaction_updated <- cellchat_db_interaction_updated %>% distinct(complex_interaction, .keep_all = TRUE)
# log_info(glue("Number of interactions after removing duplicates: {nrow(cellchat_db_interaction_updated)}"))

# # Updating complexes
# cellchat_db_complex_updated <- rbind(cellchat_db_complex, new_cellchat_db_complex)
# log_info(glue("Number of interactions: {nrow(cellchat_db_complex)}"))
# cellchat_db_complex_updated <- cellchat_db_complex_updated %>% distinct()
# log_info(glue("Number of interactions: {nrow(cellchat_db_complex_updated)}"))

# # Make sure that the complexes are unique (rowname might not be equal to combined complexes)
# cellchat_db_complex_updated <- cellchat_db_complex_updated %>%
#     unite("merged_subunits", subunit_1:subunit_5, sep = "", remove = FALSE) %>%
#     distinct(merged_subunits, .keep_all = TRUE) %>%
#     select(merged_subunits)
# log_info(glue("Number of interactions: {nrow(cellchat_db_complex_updated)}"))

# # Updating cofactors
# cellchat_db_cofactor_updated <- rbind(cellchat_db_cofactor, new_cellchat_db_cofactor)
# log_info(glue("Number of interactions: {nrow(cellchat_db_cofactor)}"))
# cellchat_db_cofactor_updated <- cellchat_db_cofactor_updated %>% distinct()
# log_info(glue("Number of interactions: {nrow(cellchat_db_cofactor_updated)}"))

# # Updating gene info
# cellchat_db_geneinfo_updated <- rbind(cellchat_db_geneinfo, new_cellchat_db_geneinfo)
# log_info(glue("Number of interactions: {nrow(cellchat_db_geneinfo)}"))
# cellchat_db_geneinfo_updated <- cellchat_db_geneinfo_updated %>% distinct()
# log_info(glue("Number of interactions: {nrow(cellchat_db_geneinfo_updated)}"))

# log_info(glue("Number of interactions originally: {nrow(cellchat_db_interaction)}"))
# log_info(glue("Number of interactions after update: {nrow(cellchat_db_interaction_updated)}"))

# log_info(glue("Number of interactions originally: {nrow(cellchat_db_complex)}"))
# log_info(glue("Number of interactions after update: {nrow(cellchat_db_complex_updated)}"))

# CellChatDB <- list()
# CellChatDB$interaction <- cellchat_db_interaction_updated
# CellChatDB$complex <- cellchat_db_complex_updated
# CellChatDB$cofactor <- cellchat_db_cofactor_updated
# CellChatDB$geneInfo <- cellchat_db_geneinfo_updated
saveRDS(cellchatDB_Omni, file = glue("{args$output_dir}/cellchat_db.rds"))
