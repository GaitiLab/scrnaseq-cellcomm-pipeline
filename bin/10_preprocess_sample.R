#!/usr/local/bin/_entrypoint.sh Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
options(Seurat.object.assay.version = "v4")
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, Seurat)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Preprocessing of individual samples",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "01_prepare_data",
        "03_preprocessed_objects",
        "seurat"
    )
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

params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    params$sample_id <- "Sample_2"
    params$input_file <- file.path(
        "output",
        "cci_pipeline",
        "01_prepare_data",
        "02_intermediate_objects",
        paste0(params$sample_id, ".rds")
    )
    params$annot <- "seurat_annotations"
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
#     params$input_file,
#     required_file_extension = "rds"
# )
# if (!checked_path) {
#     stop("Given input file is not a valid path.")
# }

# ---- Workflow ----
seurat_obj <- readRDS(params$input_file) |>
    scrnaseq.cellcomm::prepareData(
        annot = params$annot,
        min_cells = params$min_cells
    )
log_info("Loaded Seurat object & normalized data.")

GaitiLabUtils::create_dir(file.path(params$output_dir, "seurat"))
saveRDS(
    seurat_obj,
    file.path(params$output_dir, "seurat", paste0(params$sample_id, ".rds"))
)
log_info("Saved Seurat object.")

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
