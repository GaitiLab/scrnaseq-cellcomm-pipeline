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
        description = "Combine aggregation by sample (p-value combi) + aggregation by patient",
    )
    parser$add_argument("--interactions_by_patient", default = "", help = "Result of 402b aggregation by patient")
    parser$add_argument("--interactions_by_sample", default = "", help = "Result of 402a aggregation by sample")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    run_name <- "CCI_CellClass_L4_w_agg"
    args$output_dir <- glue("{here::here()}/output/{run_name}/402_aggregation")
    args$interactions_by_patient <- glue("{here::here()}/output/{run_name}/402_aggregation/402_patient_interactions_mvoted_w_filters.rds")
    args$interactions_by_sample <- glue("{here::here()}/output/{run_name}/402_aggregation/402_samples_interactions_aggregated.rds")
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
log_info("Load interactions aggregated by patient...")
interactions_mvoted <- readRDS(args$interactions_by_patient)

log_info("Load interactions aggregated by sample...")
samples_agg <- readRDS(args$interactions_by_sample)

log_info("Combine...")
combi <- merge(samples_agg, interactions_mvoted, by = c("Region_Grouped", "complex_interaction", "source_target"), all = TRUE) %>%
    distinct() %>%
    filter(!is.na(lenient_condition))

log_info("Filtering results...")
obj_filtered <- combi %>%
    filter(!is.na(lenient_condition), lenient_condition, pval < 0.05) %>%
    select(Region_Grouped, source_target, complex_interaction, lenient_condition_n_patients, lenient_condition_n_samples, lenient_condition_samples, lenient_condition_patients, pval) %>%
    distinct() %>%
    as.data.frame()

log_info("Save results...")
saveRDS(combi, glue("{args$output_dir}/402_interactions_combi_agg.rds"))
saveRDS(obj_filtered, file = glue("{args$output_dir}/402_interactions_combi_agg_filtered.rds"))


log_info("COMPLETED!")
