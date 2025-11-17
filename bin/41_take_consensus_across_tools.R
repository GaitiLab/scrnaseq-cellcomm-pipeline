#!/usr/local/bin/_entrypoint.sh Rscript

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, duckplyr)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Determine consensus for the CCI results of a sample",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "03_consensus",
        "intermediate_objects"
    )
)
parser$add_argument(
    "-a",
    "--alpha",
    type = "numeric",
    default = 0.05,
    help = "Significance threshold"
)
parser$add_argument(
    "-id",
    "--sample_id",
    type = "character",
    default = 1,
    help = "Sample id"
)

parser$add_argument(
    "--cellchat_obj",
    type = "character",
    default = "",
    help = "CellChat object"
)
parser$add_argument(
    "--liana_obj",
    type = "character",
    default = "",
    help = "LIANA object"
)
parser$add_argument(
    "--cell2cell_obj",
    type = "character",
    default = "",
    help = "Cell2Cell object"
)
parser$add_argument(
    "--cpdb_obj",
    type = "character",
    default = "",
    help = "CPDB object"
)

params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs

    params$cci_dir <- file.path(
        "output",
        "02_run_cci"
    )

    params$output_dir <- file.path(
        "output",
        "03_consensus",
        "intermediate_objects"
    )
    params$sample_id <- "BRCA2_1121251"

    params$cell2cell_obj <- file.path(
        params$cci_dir,
        "01_cell2cell",
        "02_formatted",
        paste("cell2cell", params$sample_id, "postproc.rds", sep = "__")
    )

    params$cellchat_obj <- file.path(
        params$cci_dir,
        "02_cellchat",
        "02_formatted",
        paste("cellchat", params$sample_id, "postproc.rds", sep = "__")
    )

    params$cpdb_obj <- file.path(
        params$cci_dir,
        "03_cellphonedb",
        "02_formatted",
        paste("cpdb", params$sample_id, "postproc.rds", sep = "__")
    )
    params$liana_obj <- file.path(
        params$cci_dir,
        "04_liana",
        "02_formatted",

        paste("liana", params$sample_id, "postproc.rds", sep = "__")
    )
    params$output_dir <- file.path(
        "output",
        "03_consensus",
        "intermediate_objects"
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

arg_paths <- c(
    params$cellchat_obj,
    params$liana_obj,
    params$cell2cell_obj,
    params$cpdb_obj
)

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

# ---- Load data ----
cols_common <- c(
    "source_target",
    "pval",
    "interaction_score",
    "complex_interaction",
    "method",
    "sample_id"
)

interactions_df <- arg_paths |>
    purrr::map_dfr(\(path) {
        readRDS(path) |> dplyr::select(dplyr::all_of(cols_common))
    }) |>
    duckplyr::as_duckdb_tibble()
log_info("Data loaded successfully.")

# ---- Data wrangling ----
interactions_detected_in_methods_df <- scrnaseq.cellcomm::TakeConsensus(
    interactions_df,
    alpha = params$alpha
)
log_info("Performed consensus.")

saveRDS(
    interactions_detected_in_methods_df,
    file.path(
        params$output_dir,
        paste(params$sample_id, "interactions_mvoted.rds", sep = "__")
    )
)
log_info("Output saved.")

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
