#!/usr/local/bin/_entrypoint.sh Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, duckplyr)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Post-processing LIANA results",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "02_run_cci",
        "04_liana",
        "02_formatted"
    )
)
parser$add_argument(
    "--input_interactions",
    type = "character",
    default = "",
    help = "Path to LIANA results"
)
parser$add_argument(
    "--sample_id",
    type = "character",
    default = NULL,
    help = "Sample ID"
)
parser$add_argument(
    "--ref_db",
    type = "character",
    help = "Path to interactions database",
    default = "data/interactions_db/ref_db.rds"
)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    # Provide arguments here for local runs
    cci_dir <- file.path(
        "/cluster/projects/gaitigroup/Users/Joan/breast/breast_scrnaseq/07_output/CCI/BRCA2_BSO_Neg_CellClass_L2.5/02_run_cci"
    )
    params$sample_id <- "BRCA2_1235993"
    params$ref_db <- "/cluster/projects/gaitigroup/Pipelines/scrnaseq-cellcomm-pipeline/assets/interactions_db/ref_db.rds"
    params$output_dir <- file.path("output", "testing")
    params$input_interactions <- file.path(
        cci_dir,
        "04_liana",
        "01_raw",
        paste0(paste("liana", paste0(params$sample_id, ".rds"), sep = "__"))
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
# arg_paths <- c(params$input_interactions, params$ref_db)
# checked_filepaths <- data.frame(
#     path = arg_paths,
#     required_file_extension = rep("rds", 2)
# ) |>
#     purrr::pmap_lgl(GaitiLabUtils::is_valid_path) |>
#     setNames(nm = arg_paths)
# if (!all(checked_filepaths)) {
#     stop(
#         "Not all valid paths, please check the following inputs\n",
#         paste(names(checked_filepaths)[!checked_filepaths], collapse = "\n")
#     )
# }

# ---- Load data ----
cci_obj <- readRDS(params$input_interactions)
log_info("Loaded LIANA object.")

ref_db <- readRDS(params$ref_db) |>
    dplyr::select(
        complex_interaction,
        interaction
    )
log_info("Loaded ref. database.")
# ---- Data wrangling ----
interactions_df <- scrnaseq.cellcomm::FormatLIANA(
    cci_obj = cci_obj,
    sample_id = params$sample_id,
    ref_db = ref_db
)
log_info("Standardized format of LIANA results...")

saveRDS(
    interactions_df,
    file.path(
        params$output_dir,
        paste("liana", params$sample_id, "postproc.rds", sep = "__")
    )
)
log_info("Saved LIANA results as RDS.")

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
