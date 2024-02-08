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
        description = "Get metadata",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
    args$cellchat_db <- "001_data_local/interactions_db_v2/cellchat_liana_format.rds"
    args$cpdb_db <- "001_data_local/interactions_db_v2/cpdbv5_liana_format.rds"
    args$liana_db <- "001_data_local/interactions_db_v2/liana_db_base.rds"

    args$interactions_db <- "001_data_local/interactions_db_v2/liana_db.rds"
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
# BEFORE ANY FILTERING
cellchat_db <- readRDS(args$cellchat_db) %>%
    mutate(
        method = "CellChat extracted",
        interaction = glue("{source_genesymbol} - {target_genesymbol}")
    ) %>%
    pull(interaction)
cpdb_db <- readRDS(args$cpdb_db) %>%
    mutate(method = "CellphoneDB extracted", interaction = glue("{source_genesymbol} - {target_genesymbol}")) %>%
    pull(interaction)
liana_db_base <- readRDS(args$liana_db) %>% mutate(interaction = glue("{source_genesymbol} - {target_genesymbol}"))

liana_consensus <- liana_db_base %>%
    filter(method == "LIANA consensus") %>%
    pull(interaction)

liana_rami <- liana_db_base %>%
    filter(method == "LIANA Ramilowski2015") %>%
    pull(interaction)


interactions <- list(
    cellchat_db = cellchat_db,
    cpdb_db = cpdb_db,
    liana_db_consensus = liana_consensus,
    liana_db_rami = liana_rami
)

create_venndiagram(
    x = interactions,
    filename = glue("{args$output_dir}/compare_dbs.png"),
    category.names = c("CellChat", "CellphoneDB", "LIANA consensus", "LIANA Ramilowski2015")
)

log_info(glue("Number of interactions in CellChat: {length(cellchat_db)}"))
log_info(glue("Number of interactions in CellphoneDB: {length(cpdb_db)}"))
log_info(glue("Number of interactions in LIANA consensus: {length(liana_consensus)}"))
log_info(glue("Number of interactions in LIANA Ramilowski2015: {length(liana_rami)}"))


# AFTER FILTERING
interactions_db <- readRDS(args$interactions_db)

cellchat_db <- interactions_db %>%
    filter(method == "CellChat extracted") %>%
    pull(complex_interaction)

cpdb_db <- interactions_db %>%
    filter(method == "CellphoneDB extracted") %>%
    pull(complex_interaction)

liana_db_consensus <- interactions_db %>%
    filter(method == "LIANA consensus") %>%
    pull(complex_interaction)

liana_db_rami <- interactions_db %>%
    filter(method == "LIANA Ramilowski2015") %>%
    pull(complex_interaction)

log_info(glue("Number of interactions in CellChat: {length(cellchat_db)}"))
log_info(glue("Number of interactions in CellphoneDB: {length(cpdb_db)}"))
log_info(glue("Number of interactions in LIANA consensus: {length(liana_consensus)}"))
log_info(glue("Number of interactions in LIANA Ramilowski2015: {length(liana_rami)}"))

interactions <- list(
    cellchat_db = cellchat_db,
    cpdb_db = cpdb_db,
    liana_db_consensus = liana_db_consensus,
    liana_db_rami = liana_db_rami
)

create_venndiagram(
    x = interactions,
    filename = glue("{args$output_dir}/compare_dbs_post_filtering.png"),
    category.names = c("CellChat", "CellphoneDB", "LIANA consensus", "LIANA Ramilowski2015")
)
