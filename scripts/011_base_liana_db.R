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
        description = "Get LIANA database (Ramilowski2015 + Consensus)",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/001_data_local/interactions_db_v2/")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Base Database (LIANA)
log_info("Loading LIANA database: Consensus...")
liana_consensus <- liana::select_resource("Consensus")[[1]] %>% mutate(method = "LIANA consensus")

log_info("Loading LIANA database: Ramilowski2015...")
liana_ramilowski <- liana::select_resource("Ramilowski2015")[[1]] %>% mutate(method = "LIANA Ramilowski2015")

log_info("Account for missing columns...")
missing_liana <- setdiff(colnames(liana_consensus), colnames(liana_ramilowski))
missing_ramilowski <- setdiff(colnames(liana_ramilowski), colnames(liana_consensus))

liana_consensus[missing_ramilowski] <- ""
liana_ramilowski[missing_liana] <- ""

log_info("Combining the two LIANA databases...")
liana_db <- rbind(liana_consensus, liana_ramilowski)

log_info("Saving LIANA database...")
saveRDS(liana_db, glue("{args$output_dir}/liana_db_base.rds"))

log_info("COMPLETED!")
