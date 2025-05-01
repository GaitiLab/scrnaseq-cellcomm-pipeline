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
    description = "Preprocessing of individual samples",
    default_output_dir = "output/100_preprocessing",
    default_log_file = NULL,
    default_log_dir = "output"
)
parser$add_argument(
    "--input_file",
    type = "character",
    default = NULL,
    help = "Path to Seurat object"
)
parser$add_argument(
    "--annot",
    type = "character",
    default = "CellClass_L1",
    help = "Annotation to use for filtering"
)
parser$add_argument(
    "-n",
    "--min_cells",
    type = "integer",
    default = 5,
    help = "Minimum number of cells required in each cell group for cell-cell communication (default=5)"
)
parser$add_argument(
    "--sample_id",
    type = "character",
    default = NULL,
    help = "Sample ID"
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
options(Seurat.object.assay.version = "v4")

log_info("Prepare data...")
scrnaseq.cellcomm::prepare_data(
    input_file = params$input_file,
    annot = params$annot,
    output_dir = params$output_dir,
    min_cells = params$min_cells,
    sample_id = params$sample_id
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
