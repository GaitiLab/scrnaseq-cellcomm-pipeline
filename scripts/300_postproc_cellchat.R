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
        description = "Post-processing CellChat",
    )
    parser$add_argument("--input_interactions", type = "character", default = NULL, help = "Directory with CellChat results")
    parser$add_argument("--ref_db", type = "character", default = NULL, help = "Path to interactions database")
    parser$add_argument("--sample_id", type = "character", default = NULL, help = "Sample ID")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/300_postproc_cellchat")
    args$input_interactions <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/200_cci_cellchat/cellchat__6514_enhancing_border.rds")
    args$ref_db <- glue("{here::here()}/001_data_local/interactions_db_v2/ref_db.rds")
    args$sample_id <- "6514_enhancing_border"
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
pacman::p_load(reshape2, tidyverse)

log_info("Load reference database...")
ref_db <- readRDS(args$ref_db) %>%
    select(
        simple_interaction,
        complex_interaction, interaction
    )

log_info("Load data...")
interactions <- readRDS(args$input_interactions) %>%
    left_join(ref_db, by = "interaction") %>%
    select(-interaction) %>%
    unite(source_target, source, target, sep = "__") %>%
    rename(proba_score = proba) %>%
    mutate(method = "CellChatv2", Sample = args$sample_id)

log_info("Save output...")
saveRDS(interactions, glue("{args$output_dir}/cellchat__{args$sample_id}__postproc.rds"))
log_info("COMPLETED!")
