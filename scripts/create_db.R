# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Create a database",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/data/interactions_db_v2_2")
    args$source_cpdb_dir <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/scrnaseq-cellcomm-pipeline/000_misc_local/references/cellphonedb_v5.0.0"
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
scrnaseq.cellcomm::create_db(
    output_dir = args$output_dir,
    source_cpdb_dir = args$source_cpdb_dir
)
