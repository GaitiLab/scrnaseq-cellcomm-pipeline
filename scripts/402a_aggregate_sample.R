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
        description = "Aggregate interactions by sample (combine p-values)",
    )
    parser$add_argument("--input_file", type = "character", help = "Input from 401_combine_samples.R, name should be 401_samples_interactions_agg_rank.rds", default = "")
    parser$add_argument("--condition_varname", default = "Region_Grouped", type = "character", help = "Grouping variable")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CCI_CellClass_L1_w_agg/402_aggregation")
    args$input_file <- glue("{here::here()}/output/CCI_CellClass_L1_w_agg/401_combine_samples/401_samples_interactions_agg_rank.rds")
    args$condition_varname <- "Region_Grouped"
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
log_info("Load additional libraries...")
pacman::p_load(survcomp)

log_info("Load interactions and combine p-values...")
obj <- readRDS(args$input_file) %>%
    ungroup() %>%
    select(!!sym(args$condition_varname), complex_interaction, pval, source_target) %>%
    group_by(!!sym(args$condition_varname), source_target, complex_interaction) %>%
    # Combine p-values across samples
    mutate(pval = combine.test(pval))

log_info("Save results...")
saveRDS(obj, glue("{args$output_dir}/402_samples_interactions_aggregated.rds"))

log_info("COMPLETED!")
