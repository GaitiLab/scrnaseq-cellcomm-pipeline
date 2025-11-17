#!/usr/local/bin/_entrypoint.sh Rscript

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, duckplyr)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Post-process Cell2Cell results",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "02_run_cci",
        "01_cell2cell",
        "02_formatted"
    )
)
parser$add_argument(
    "-is",
    "--input_interactions_scores",
    type = "character",
    default = NULL,
    help = "CSV file with interaction scores from cell2cell"
)
parser$add_argument(
    "-ip",
    "--input_interactions_pval",
    type = "character",
    default = NULL,
    help = "CSV file with pvalues from cell2cell"
)
parser$add_argument(
    "-s",
    "--sample_id",
    type = "character",
    default = NULL,
    help = "Sample ID"
)
parser$add_argument(
    "--ref_db",
    type = "character",
    default = NULL,
    help = "Path to interactions database"
)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs

    cci_dir <- file.path(
        "/cluster/projects/gaitigroup/Users/Joan/breast/breast_scrnaseq/07_output/CCI/BRCA2_BSO_Neg_CellClass_L2/02_run_cci"
    )
    params$sample_id <- "BRCA2_1121251"
    params$ref_db <- "/cluster/projects/gaitigroup/Pipelines/scrnaseq-cellcomm-pipeline/assets/interactions_db/ref_db.rds"
    params$output_dir <- file.path("output_local")
    params$input_interactions_scores <- file.path(
        cci_dir,
        "01_cell2cell",
        "01_raw",
        paste(
            "cell2cell",
            params$sample_id,
            "interaction_scores.csv",
            sep = "__"
        )
    )
    params$input_interactions_pval <- file.path(
        cci_dir,
        "01_cell2cell",
        "01_raw",
        paste(
            "cell2cell",
            params$sample_id,
            "pvalues.csv",
            sep = "__"
        )
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
arg_paths <- c(
    params$input_interactions_pval,
    params$input_interactions_scores,
    params$ref_db
)
checked_filepaths <- data.frame(
    path = arg_paths,
    required_file_extension = c(rep("csv", 2), "rds")
) |>
    purrr::pmap_lgl(GaitiLabUtils::is_valid_path) |>
    setNames(nm = arg_paths)
if (!all(checked_filepaths)) {
    stop(
        "Not all valid paths, please check the following inputs\n",
        paste(names(checked_filepaths)[!checked_filepaths], collapse = "\n")
    )
}

# ---- Load data ----
log_info("Load data...")
c2c_pvalues_df <- data.table::fread(params$input_interactions_pval)
c2c_interaction_scores_df <- data.table::fread(params$input_interactions_scores)

ref_db <- readRDS(params$ref_db) |>
    dplyr::select(
        complex_interaction,
        interaction
    )

# ---- Data wrangling ----
log_info("Standardize format of Cell2Cell results...")
interactions_combined_df <- scrnaseq.cellcomm::FormatCell2CellWrapper(
    input_interactions_scores = c2c_interaction_scores_df,
    input_interactions_pval = c2c_pvalues_df,
    sample_id = params$sample_id,
    ref_db = ref_db
)

log_info("Save outputs...")
saveRDS(
    interactions_combined_df,
    file = file.path(
        params$output_dir,
        paste("cell2cell", params$sample_id, "postproc.rds", sep = "__")
    )
)

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
