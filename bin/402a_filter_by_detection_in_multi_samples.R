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
    description = "Post-filtering/formatting",
    default_output_dir = "output",
    default_log_file = NULL,
    default_log_dir = "output"
)
parser$add_argument(
    "--input_file",
    type = "character",
    default = NULL,
    help = "Path to '401_samples_interactions_mvoted.rds' file"
)
parser$add_argument(
    "--min_patients",
    type = "integer",
    default = 2,
    help = "Minimum number of patients for an interaction to be kept"
)
parser$add_argument(
    "--condition_var",
    type = "character",
    help = "Name of condition variable",
    default = "Condition_dummy"
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
logr <- init_logging(
    log_level = params$log_level,
    log_file = NULL
)
obj_logger <- init_obj_logging(
    log_file = NULL
)
log_info(ifelse(
    interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Parameters:")
log_object(params_ls_to_df(params))

scrnaseq.cellcomm::filter_by_detection_in_multi_samples(
    input_file = params$input_file,
    min_patients = params$min_patients,
    condition_var = params$condition_var,
    output_dir = params$output_dir
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
