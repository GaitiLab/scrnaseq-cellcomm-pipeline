#!/usr/local/bin/_entrypoint.sh Rscript

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, duckplyr)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Aggregate interactions by sample (combine p-values)",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "04_aggregation"
    )
)
parser$add_argument(
    "--input_file",
    type = "character",
    help = "Input file samples_interactions_agg_rank.rds",
    default = ""
)
parser$add_argument(
    "--method",
    type = "character",
    help = "Correction method",
    default = 'BH'
)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs

    # Provide arguments here for local runs
    params$input_file <- file.path(
        "output",
        "03_consensus",
        "samples_interactions_agg_rank.rds"
    )
}

create_dir(params$output_dir)
# Set up logging
log_file <- set_logfile(params)
logr <- init_logging(log_level = params$log_level, log_file = log_file)
obj_logger <- init_obj_logging(log_file = log_file)
log_info(ifelse(
    interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Parameters:")
log_object(params_ls_to_df(params))

# ---- Check arguments ----
checked_path <- is_valid_path(
    params$input_file,
    required_file_extension = "rds"
)
if (!checked_path) {
    stop("Given input file is not a valid path.")
}

# ---- Workflow ----
obj <- readRDS(params$input_file)
log_info("Loaded dataframe.")

obj_agg <- obj |>
    scrnaseq.cellcomm::AggregateSamples(method = params$method)
log_info("Aggregated samples.")

saveRDS(
    obj_agg,
    file.path(params$output_dir, "aggregation_samples.rds")
)
log_info("Saved dataframe.")

log_info("Finished")

log_info("Session Info")
log_object(sessionInfo())

if (!is.null(params$nf_process_id)) {
    write_versions_yml(
        c("scrnaseq.cellcomm", pacman::p_loaded()),
        task_id = params$nf_process_id,
        outdir = params$output_dir
    )
}
