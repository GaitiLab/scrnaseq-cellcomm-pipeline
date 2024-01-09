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
    args$output_dir <- glue("/Users/joankant/Library/CloudStorage/OneDrive-SharedLibraries-UHN/Wu, Yiyan - Spatial_GBM/Analysis/CCI")
    args$sheet_name <- "CCI_CellClass_L2"
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
require(xlsx)

log_info("Load data...")
obj <- readRDS(glue("{here::here()}/output/{args$sheet_name}/400_consensus/400_samples_interactions_mvoted_w_filters.rds"))

obj_filtered <- obj %>%
    filter(lenient_region_pair) %>%
    select(Region_Grouped, source_target, source, target, complex_interaction, lenient_voting_samples) %>%
    distinct() %>%
    as.data.frame()

filename <- glue("{args$output_dir}/interactions.xlsx")
write.xlsx(obj_filtered, filename,
    sheetName = args$sheet_name,
    col.names = TRUE, row.names = FALSE, append = TRUE
)
