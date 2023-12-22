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
        description = "Test",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
# create_dir(args$output_dir)

# Load additional libraries

liana_db <- readRDS("001_data_local/interactions_db_v2/liana_db__20231211.rds")
write.csv(liana_db, "001_data_local/interactions_db_v2/liana_db__20231211.csv", row.names = FALSE)

new <- read.csv("001_data_local/interactions_db_v2/liana_db__20231211.csv")

old <- read.csv("001_data_local/interactions_db/custom_liana.csv")

liana_old <- readRDS("001_data_local/interactions_db/custom_liana.rds")
