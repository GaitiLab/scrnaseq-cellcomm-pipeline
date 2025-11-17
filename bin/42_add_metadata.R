#!/usr/local/bin/_entrypoint.sh Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, duckplyr)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Combine samples",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "03_consensus"
    )
)
parser$add_argument(
    "--input_dir",
    help = "Path to input directory",
    type = "character"
)
parser$add_argument(
    "--meta_df",
    type = "character",
    help = "Path to meta_df (RDS file)",
    default = NULL
)
parser$add_argument(
    "--suffix",
    type = "character",
    default = "interactions_mvoted",
    help = "Pattern to look for in `input_dir`, i.e. 'interactions_mvoted' or 'interactions_agg_rank'"
)
parser$add_argument(
    "--sample_var",
    type = "character",
    help = "Name of sample variable (default = 'sample_id')",
    default = "sample_id"
)
parser$add_argument(
    "--condition_var",
    type = "character",
    help = "Name of condition variable",
    default = "Condition_dummy"
)
parser$add_argument(
    "--patient_var",
    type = "character",
    help = "Name of patient variable",
    default = "sample_id"
)

params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    params$sample_var <- "ID"
    params$condition_var <- "disease_status"
    params$patient_var <- "patient"
    params$input_dir <- file.path(
        "output",
        "PR_GBM",
        "03_consensus",
        "intermediate_objects"
    )
    # params$meta_df <- "/cluster/projects/gaitigroup/Users/Jiaoyi/breast_scrnaseq/07_output/CCI/BRCA2_BSO_Neg_CellClass_L2/01_prepare_data/01_misc/BRCA2_neg_BSO_res05_scVI_w_anno__metadata.rds"
    params$meta_df <- "/cluster/projects/gaitigroup/Users/Yiyan/Reanalysis/07_output/PR_reanalysis/CCI/01_prepare_data/01_misc/combined_matched_PR1__metadata.rds"
    params$suffix <- "interactions_mvoted"
}

create_dir(params$output_dir)
# Set up logging
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
arg_paths <- c(params$input_dir, params$meta_df)

checked_filepaths <- data.frame(
    path = arg_paths,
    expected_type = c("dir", "file"),
    required_file_extension = c(NA, "rds")
) |>
    purrr::pmap_lgl(GaitiLabUtils::is_valid_path) |>
    setNames(nm = as.character(arg_paths))
if (!all(checked_filepaths)) {
    stop(
        "Not all valid paths, please check the following inputs\n",
        paste(names(checked_filepaths)[!checked_filepaths], collapse = "\n")
    )
}
# ---- Workflow ----
paths <- list.files(
    params$input_dir,
    pattern = paste0(params$suffix, ".rds"),
    full.names = TRUE
)
df <- paths |>
    purrr::map_dfr(readRDS)
log_info(
    paste(
        "Loaded all dataframes",
        paste0("(n=", length(paths), ")"),
        "with filename ending with",
        paste0(params$suffix, ".rds")
    )
)
meta_df <- readRDS(params$meta_df)
log_info("Loaded metadata.")

meta_clean_df <- meta_df |>
    scrnaseq.cellcomm::FormatMetadata(
        patient_var = params$patient_var,
        condition_var = params$condition_var,
        sample_var = params$sample_var
    )


df <- df |> left_join(meta_clean_df)
log_info("Combined all samples and added metadata.")

saveRDS(
    df,
    file = file.path(
        params$output_dir,
        paste("samples", paste0(params$suffix, ".rds"), sep = "_")
    )
)
log_info("Saved output as rds.")

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
