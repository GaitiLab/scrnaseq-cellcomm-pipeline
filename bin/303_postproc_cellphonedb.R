#!/usr/local/bin/_entrypoint.sh Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)


# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
# Comment/uncomment depending on whether you have an internal package based on the 'R' directory created with usethis::create_package()
# devtools::load_all("./", export_all = FALSE)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Post-processing of CellPhoneDB results",
    default_output_dir = "output/303_postproc_cpdb",
    default_log_file = NULL,
    default_log_dir = "output"
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
parser$add_argument(
    "--task_id",
    type = "character",
    help = "Task ID from Nextflow, only needed in Nextflow pipeline",
    default = NULL
)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
}

create_dir(params$output_dir)
# Set up logging
logr <- init_logging(log_level = params$log_level, log_file = NULL)
obj_logger <- init_obj_logging(log_file = NULL)
log_info(ifelse(
    interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Parameters:")
log_object(params_ls_to_df(params))

log_info("Standardize format of CPDB results...")
scrnaseq.cellcomm::format_cpdb(
    interaction_scores = params$interaction_scores,
    pval = params$pval,
    sign_means = params$sign_means,
    means = params$means,
    output_dir = params$output_dir,
    sample_id = params$sample_id,
    ref_db = params$ref_db
)

log_info("Finished")

log_info("Session Info")
log_object(sessionInfo())

if (!is.null(params$task_id)) {
    write_versions_yml(
        c("scrnaseq.cellcomm", pacman::p_loaded()),
        task_id = params$task_id,
        outdir = params$output_dir
    )
}
