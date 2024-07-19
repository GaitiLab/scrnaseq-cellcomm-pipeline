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
        description = "Create ExcelSheet with interactions from post-filtered data"
    )
    parser$add_argument("--interactions_agg_integration", type = "character", help = "path to 402c_filtering_aggregated_res.rds")
    parser$add_argument("--condition_var", type = "character", help = "Condition variablem e.g. mutation or region")
    parser$add_argument("--alpha", type = "numeric", help = "Alpha for additional filtering (default = 1.01, e.g. keeping all)", default = 1.01)
    parser$add_argument("--output_name", type = "character", default = "interactions_summary", help = "filename without extension for saving interactions in an Excel file.")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- "output/test_individual_scripts/"
    args$interactions_agg_integration <- "output/test_individual_scripts/402_aggregation_and_filtering/402c_filtering_aggregated_res.rds"
    args$output_name <- "interactions_summary"
    args$condition_var <- "Condition"
    args$alpha <- 0.05

    args$interactions_agg_integration <- "output/LP_IMM_perSample/402_aggregation_and_filtering/402c_filtering_aggregated_res.rds"
    args$condition_var <- "Mutation"
    args$output_dir <- "output/LP_IMM_perSample"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

log_info("Format and save as xlsx...")
scrnaseq.cellcomm::save_as_xlsx(
    interactions_agg_integration = args$interactions_agg_integration,
    condition_var = args$condition_var,
    alpha = args$alpha,
    output_dir = args$output_dir,
    output_name = args$output_name
)
log_info("Finished!")
