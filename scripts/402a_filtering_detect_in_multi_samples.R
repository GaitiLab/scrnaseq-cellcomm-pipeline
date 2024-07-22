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
        description = "Post-filtering/formatting",
        default_output = "output/402_aggregation_and_filtering"
    )
    parser$add_argument("--input_file",
        type = "character", default = NULL,
        help = "Path to '401_samples_interactions_mvoted.rds' file"
    )
    parser$add_argument("--min_patients", type = "integer", default = 2, help = "Minimum number of patients for an interaction to be kept")
    parser$add_argument("--condition_var", type = "character", help = "Name of condition variable", default = "Condition_dummy")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- "output/test_individual_scripts/402_aggregation_and_filtering"
    args$input_file <- "output/test_individual_scripts/401_combine_samples/401_samples_interactions_mvoted.rds"
    args$condition_var <- "Condition"
    args$min_patients <- 1
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
scrnaseq.cellcomm::filter_by_detection_in_multi_samples(
    input_file = args$input_file,
    min_patients = args$min_patients,
    condition_var = args$condition_var,
    output_dir = args$output_dir
)
