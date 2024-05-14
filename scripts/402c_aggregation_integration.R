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
        description = "Combine aggregation by sample (p-value combi) + aggregation by patient",
    )
    parser$add_argument("--input_dir", default = "", help = "Directory called 402_aggregation")
    parser$add_argument("--interactions_agg_binarized", default = "", help = "Path to 402a_aggregation_binarized.rds")
    parser$add_argument("--interactions_agg_continuous", default = "", help = "Path to 402b_aggregation_continuous.rds")
    parser$add_argument("--condition_varname", default = "", help = "Condition, e.g. mutation, region")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    run_name <- "CCI_CellClass_L2_2_reassigned_samples_confident_only"
    args$output_dir <- glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/GBM_CCI_Analysis/output/{run_name}/402_aggregation")
    args$input_dir <- glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/GBM_CCI_Analysis/output/{run_name}/402_aggregation")
    args$condition_varname <- "Region"
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
if (file.exists(args$input_dir)) {
    log_info("Load interactions after aggregation: binarized...")
    interactions_binarized <- readRDS(glue("{args$input_dir}/402a_aggregation_binarized.rds"))

    log_info("Load interactions after aggregation: continuous...")
    interactions_continuous <- readRDS(glue("{args$input_dir}/402b_aggregation_continuous.rds"))
} else {
    log_info("Load interactions after aggregation: binarized...")
    interactions_binarized <- readRDS(args$interactions_agg_binarized)

    log_info("Load interactions after aggregation: continuous...")
    interactions_continuous <- readRDS(args$interactions_agg_continuous)
}

log_info("Combine...")
combi <- merge(interactions_continuous, interactions_binarized, by = c(args$condition_varname, "complex_interaction", "source_target"), all = TRUE) %>%
    distinct() %>%
    filter(!is.na(lenient_condition))

log_info("Save results...")
saveRDS(combi, glue("{args$output_dir}/402c_aggregation_integration.rds"))

log_info("COMPLETED!")
