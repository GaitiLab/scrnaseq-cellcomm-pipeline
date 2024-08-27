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
        description = "Aggregate interactions by sample (combine p-values)",
        default_output = "output/402_aggregation_and_filtering"
    )
    parser$add_argument("--input_file", type = "character", help = "Input from file 401_samples_interactions_agg_rank.rds, name should be 401_samples_interactions_agg_rank.rds", default = "")
    parser$add_argument("--condition_var", default = "Condition_dummy", type = "character", help = "Grouping variable")
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

log_info("Aggregate all samples with ranked interactions...")
scrnaseq.cellcomm::aggregate_samples(
    input_file = args$input_file,
    output_dir = args$output_dir,
    condition_var = args$condition_var
)
log_info("Finished!")
