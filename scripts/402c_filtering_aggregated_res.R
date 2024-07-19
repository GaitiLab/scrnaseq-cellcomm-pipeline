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
        description = "Combine aggregation by sample (p-value combi) + aggregation by patient",
        default_output = "output/402_aggregation_and_filtering"
    )
    parser$add_argument("--input_dir", default = "", help = "path to directory '402_aggregation_and_filtering'")
    parser$add_argument("--interactions_agg_binarized", default = "", help = "Path to 402a_filtering_detect_in_multi_samples.rds")
    parser$add_argument("--interactions_agg_continuous", default = "", help = "Path to 402b_aggregation_samples.rds")
    parser$add_argument("--condition_var", default = "Condition_dummy", help = "Condition, e.g. mutation, region")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- "output/test_individual_scripts/402_aggregation_and_filtering"
    args$input_dir <- "output/test_individual_scripts/402_aggregation_and_filtering"
    args$interactions_agg_binarized <- ""
    args$interactions_agg_continuous <- ""
    args$condition_var <- "Condition"


    args$output_dir <- "output/LP_IMM_perSample/402_aggregation_and_filtering"
    args$input_dir <- args$output_dir
    args$condition_var <- "Mutation"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

log_info("Aggregate results w/ ranked interactions from all samples...")
scrnaseq.cellcomm::filter_agg_res(
    condition_var = args$condition_var,
    input_dir = args$input_dir,
    interactions_agg_binarized = args$interactions_agg_binarized,
    interactions_agg_continuous = args$interactions_agg_continuous,
    output_dir = args$output_dir
)
log_info("Finished!")
