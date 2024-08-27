# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Post-processing LIANA results", default_output = "output/301_postproc_liana"
    )
    parser$add_argument("--input_interactions", type = "character", default = "", help = "Path to LIANA results")
    parser$add_argument("--sample_id", type = "character", default = NULL, help = "Sample ID")
    parser$add_argument("--ref_db", type = "character", help = "Path to interactions database", default = "data/interactions_db/ref_db.rds")
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
create_dir(args$output_dir)

log_info("Standardize format of LIANA results...")
scrnaseq.cellcomm::format_liana(
    input_interactions = args$input_interactions,
    output_dir = args$output_dir,
    sample_id = args$sample_id,
    ref_db = args$ref_db
)
log_info("Finished!")
