#!/usr/local/bin/_entrypoint.sh Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, duckplyr)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Combine aggregation by sample (p-value combi) + aggregation by patient",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "04_aggregation"
    )
)

parser$add_argument(
    "--interactions_mvoted",
    default = "",
    help = "Path to 402a_filtering_detect_in_multi_samples.rds"
)
parser$add_argument(
    "--interactions_ranked",
    default = "",
    help = "Path to 402b_aggregation_samples.rds"
)

params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs

    params$interactions_ranked <- file.path(
        "output",
        "04_aggregation",
        "aggregation_samples.rds"
    )
    params$interactions_mvoted <- file.path(
        "output",
        "04_aggregation",
        "filtering_detect_in_multi_samples.rds"
    )
    params$output_dir <- file.path("output", "04_aggregation")
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
arg_paths <- c(params$interactions_mvoted, params$interactions_ranked)
checked_filepaths <- data.frame(
    path = arg_paths,
    required_file_extension = c(rep("rds", length(arg_paths)))
) |>
    purrr::pmap_lgl(GaitiLabUtils::is_valid_path) |>
    setNames(nm = arg_paths)
if (!all(checked_filepaths)) {
    stop(
        "Not all valid paths, please check the following inputs\n",
        paste(names(checked_filepaths)[!checked_filepaths], collapse = "\n")
    )
}

# ---- Workflow ----
list_of_dfs <- arg_paths |>
    purrr::map(\(x) {
        readRDS(x) |> duckplyr::as_duckdb_tibble()
    })
log_info("Loaded dataframes.")

df_combined <- list_of_dfs |> purrr::reduce(left_join)
log_info("Added pvalues + scores to mvoted dataframe.")

saveRDS(
    df_combined,
    file.path(params$output_dir, "filtering_aggregated_res.rds")
)
log_info("Saved dataframe as rds.")

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
