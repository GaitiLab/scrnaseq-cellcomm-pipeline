# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# require(GBMutils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)
if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Format CellChat DB for LIANA DB",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/001_data_local/interactions_db_v2")
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

log_info("Load CellChat database...")
cellchat_db <- CellChatDB.human
cellchat_db_interaction <- cellchat_db$interaction
rownames(cellchat_db_interaction) <- NULL

log_info("Adapt for LIANA DB...")
interactions <- cellchat_db_interaction %>%
    select(interaction_name, ligand, receptor, ligand.symbol, receptor.symbol) %>%
    mutate(
        source_genesymbol = str_replace_all(ligand.symbol, ", ", "_"),
        target_genesymbol = str_replace_all(receptor.symbol, ", ", "_")
    ) %>%
    mutate(
        source_genesymbol = ifelse(source_genesymbol == "", toupper(ligand), source_genesymbol),
        target_genesymbol = ifelse(target_genesymbol == "", toupper(receptor), target_genesymbol)
    ) %>%
    select(source_genesymbol, target_genesymbol)

log_info("Save CellChat in LIANA db format...")
saveRDS(interactions, glue("{args$output_dir}/cellchat_liana_format.rds"))

log_info("COMPLETED!")
