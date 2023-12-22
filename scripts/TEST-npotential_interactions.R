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
    args$genes <- glue("{here::here()}/output/captured_genes.rds")
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
genes <- readRDS(args$genes)

# V1 interactions DB
args$interactions_db <- glue("{here::here()}/001_data_local/interactions_db/interactions_ref.rds")
interactions_db <- readRDS(args$interactions_db) %>% filter(!is.na(genename_a) & !is.na(genename_b) & (genename_a != "") & (genename_b != ""))
log_info(glue("Number of interactions: {nrow(interactions_db)}"))
n_possible <- interactions_db %>% filter((genename_a %in% genes) & (genename_b %in% genes))
log_info(glue("Number of interactions: {nrow(n_possible)}"))
log_info(glue("% covered: {round(nrow(n_possible)/nrow(interactions_db) * 100,2)}"))

# V2 interactions DB
args$interactions_db <- "001_data_local/interactions_db_v2/liana_db_final__20231116.rds"
interactions_db <- readRDS(args$interactions_db) %>% filter(!is.na(source_genesymbol) & !is.na(target_genesymbol) & (source_genesymbol != "") & (target_genesymbol != ""))
log_info(glue("Number of interactions: {nrow(interactions_db)}"))
n_possible <- interactions_db %>% filter((source_genesymbol %in% genes) & (target_genesymbol %in% genes))

log_info(glue("Number of interactions: {nrow(n_possible)}"))
log_info(glue("% covered: {round(nrow(n_possible)/nrow(interactions_db) * 100,2)}"))
