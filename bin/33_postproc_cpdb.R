#!/usr/local/bin/_entrypoint.sh Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, duckplyr)
# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Post-processing of CellPhoneDB results",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "02_run_cci",
        "03_cellphonedb",
        "02_formatted"
    )
)
parser$add_argument(
    "--sample_id",
    default = "",
    type = "character",
    help = "Sample ID"
)
parser$add_argument(
    "--interaction_scores",
    default = "",
    type = "character",
    help = "Path to CellPhoneDB interaction scores"
)
parser$add_argument(
    "--pval",
    default = "",
    type = "character",
    help = "Path to CellPhoneDB p-values"
)
parser$add_argument(
    "--sign_means",
    default = "",
    type = "character",
    help = "Path to CellPhoneDB significant means"
)
parser$add_argument(
    "--means",
    default = "",
    type = "character",
    help = "Path to CellPhoneDB means"
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
    params$interaction_scores <- file.path(
        cci_dir,
        "03_cellphonedb",
        "01_raw",
        paste0(
            paste(
                "statistical_analysis_interaction_scores",
                params$sample_id,
                sep = "__"
            ),
            ".txt"
        )
    )

    params$pval <- file.path(
        cci_dir,
        "03_cellphonedb",
        "01_raw",
        paste0(
            paste("statistical_analysis_pvalues", params$sample_id, sep = "__"),
            ".txt"
        )
    )
    params$sign_means <- file.path(
        cci_dir,
        "03_cellphonedb",
        "01_raw",
        paste0(
            paste(
                "statistical_analysis_significant_means",
                params$sample_id,
                sep = "__"
            ),
            ".txt"
        )
    )
    params$means <- file.path(
        cci_dir,
        "03_cellphonedb",
        "01_raw",
        paste0(
            paste("statistical_analysis_means", params$sample_id, sep = "__"),
            ".txt"
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
cpdb_paths <- setNames(
    c(
        params$interaction_scores,
        params$pval,
        params$sign_means,
        params$means
    ),
    c("interaction_score", "pval", "sign_mean", "mean")
)

arg_paths <- c(cpdb_paths, params$ref_db)
checked_filepaths <- data.frame(
    path = arg_paths,
    required_file_extension = c(rep("txt", 4), "rds")
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
cpdb_dfs <- cpdb_paths |>
    map(data.table::fread, sep = "\t", header = TRUE, check.names = FALSE) |>
    setNames(nm = names(cpdb_paths))

ref_db <- readRDS(params$ref_db) |>
    dplyr::select(complex_interaction, interaction)

# ---- Data wrangling ----
log_info("Standardize format of CPDB results...")
interactions_df <- scrnaseq.cellcomm::FormatCPDBWrapper(
    cpdb_dfs = cpdb_dfs,
    sample_id = params$sample_id,
    ref_db = ref_db
)

log_info("Save output...")
saveRDS(
    interactions_df,
    file.path(
        params$output_dir,
        paste("cpdb", params$sample_id, "postproc.rds", sep = "__")
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
