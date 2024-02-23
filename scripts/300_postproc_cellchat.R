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
    args$output_dir <- glue("{here::here()}/output/test_downsampling_implementation/300_postproc_cellchat")
    args$input_interactions <- glue("{here::here()}/output/CCI_CellClass_L1_conf_malign/200_cci_cellchat/cellchat__6234_2895153_A.rds")
    args$ref_db <- glue("{here::here()}/data/interactions_db/ref_db.rds")
    args$sample_id <- "6234_2895153_A__run__1"
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
pacman::p_load(reshape2)

log_info("Load reference database...")
ref_db <- readRDS(args$ref_db) %>%
    select(
        simple_interaction,
        complex_interaction, interaction
    )

# Setting sample + run id depending on downsampling or not
split_sample_id <- str_split(args$sample_id, "__", simplify = TRUE)
sample_id <- split_sample_id[1]
run_id <- NA
if (length(split_sample_id) > 2) {
    run_id <- split_sample_id[3]
}

log_info("Load data...")
interactions <- readRDS(args$input_interactions) %>%
    dplyr::rename(interaction_score = proba) %>%
    left_join(ref_db, by = "interaction") %>%
    select(-interaction) %>%
    unite(source_target, source, target, sep = "__") %>%
    mutate(method = "CellChatv2", Sample = sample_id, run_id = run_id)

log_info("Save output...")
saveRDS(interactions, glue("{args$output_dir}/cellchat__{args$sample_id}__postproc.rds"))
log_info("COMPLETED!")
