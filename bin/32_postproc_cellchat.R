#!/usr/local/bin/_entrypoint.sh Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, duckplyr)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Post-processing CellChat",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "02_run_cci",
        "02_cellchat",
        "02_formatted"
    )
)
parser$add_argument(
    "--input_interactions",
    type = "character",
    default = "",
    help = "Directory with CellChat results"
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
    default = "data/interactions_db/ref_db.rds",
    help = "Path to interactions database (default = 'data/interactions_db/ref_db.rds')"
)
parser$add_argument(
    "--n_cores",
    type = "numeric",
    default = 1,
    help = "No. of cores to use"
)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    cci_dir <- file.path(
        "/cluster/projects/gaitigroup/Users/Joan/breast/breast_scrnaseq/07_output/CCI/BRCA2_BSO_Neg_CellClass_L2/02_run_cci"
    )
    params$sample_id <- "BRCA2_1121251"
    params$ref_db <- "/cluster/projects/gaitigroup/Pipelines/scrnaseq-cellcomm-pipeline/assets/interactions_db/ref_db.rds"
    params$output_dir <- file.path("output", "testing")
    params$input_interactions <- file.path(
        cci_dir,
        "02_cellchat",
        "01_raw",
        paste0(paste("cellchat", params$sample_id, "raw_obj.rds", sep = "__"))
    )
    sample_id <- params$sample_id
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
log_info("Load data...")
cellchat_obj <- readRDS(params$input_interactions)

ref_db <- readRDS(params$ref_db) |>
    dplyr::select(complex_interaction, interaction)

# ---- Data wrangling ----
log_info("Standardize format of CellChat results...")
cc_df <- cellchat_obj |>
    scrnaseq.cellcomm::FormatCellChat(
        sample_id = params$sample_id,
        ref_db = ref_db,
        n_cores = params$n_cores
    )

log_info("Save output...")
saveRDS(
    cc_df,
    file.path(
        params$output_dir,
        paste("cellchat", params$sample_id, "postproc.rds", sep = "__")
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
