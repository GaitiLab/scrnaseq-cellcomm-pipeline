#!/usr/local/bin/_entrypoint.sh Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Create ExcelSheet with interactions from post-filtered data",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "04_aggregation"
    )
)
parser$add_argument(
    "--interactions_agg_integration",
    type = "character",
    help = "path to filtering_aggregated_res.rds"
)
parser$add_argument(
    "--alpha",
    type = "numeric",
    help = "Alpha for additional filtering (default = 1.01, e.g. keeping all)",
    default = 1.01
)
parser$add_argument(
    "--output_name",
    type = "character",
    default = "interactions_summary",
    help = "filename without extension for saving interactions in an Excel file."
)
parser$add_argument(
    "--is_stringent",
    action = "store_true",
    default = FALSE,
    help = "Indicate whether only interactions passing the filter criteria using the stringent voting method for detection in multiple tools should be kept."
)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    params$interactions_agg_integration <- file.path(
        "output",
        "04_aggregation",
        "filtering_aggregated_res.rds"
    )
    params$alpha <- 0.05
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
# checked_path <- is_valid_path(
#     params$interactions_agg_integration,
#     required_file_extension = "rds"
# )
# if (!checked_path) {
#     stop("Given input file is not a valid path.")
# }

# ---- Workflow ----

df <- readRDS(params$interactions_agg_integration) |>
    duckplyr::as_duckdb_tibble()
log_info("Loaded dataframe.")

output_filename <- file.path(
    params$output_dir,
    paste0(params$output_name, ".xlsx")
)

if (fs::file_exists(output_filename)) {
    file.remove(output_filename)
    log_info("Removed already existing output file path.")
}

if (params$is_stringent) {
    df_filtered <- df |>
        dplyr::filter(
            pval_adj < params$alpha,
            stringent_voting_detected_in_enough_patients
        ) |>
        dplyr::distinct() |>
        as.data.frame()
    log_info(
        "Removed all interactions that did not pass using the 'stringent' voting strategy and are not significant based on RRA."
    )
} else {
    df_filtered <- df |>
        dplyr::filter(pval_adj < params$alpha) |>
        dplyr::distinct() |>
        as.data.frame()
    log_info("Removed not significant interactions based on RRA.")
}
openxlsx::write.xlsx(
    df_filtered |> as.data.frame(),
    output_filename
)
log_info("Saved dataframe as Excel file.")

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
