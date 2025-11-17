#!/usr/local/bin/_entrypoint.sh Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
options(Seurat.object.assay.version = "v4")
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, Seurat)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Remove redundant assays from Seurat object",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "01_prepare_data",
        "02_intermediate_objects"
    )
)
parser$add_argument(
    "-i",
    "--input_file",
    type = "character",
    default = NULL,
    help = "Path to Seurat object"
)

params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    params$input_file <- "/cluster/projects/gaitigroup/Users/Joan/nextflow/scrnaseq-cellcomm-pipeline/test_data/example_data.rds"
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

checked_path <- is_valid_path(
    params$input_file,
    required_file_extension = "rds"
)
if (!checked_path) {
    stop("Given input file is not a valid path.")
}

seurat_obj <- readRDS(params$input_file) |>
    scrnaseq.cellcomm::RemoveUnusedAssays()
log_info("Loaded the Seurat object & reduced its size.")

saveRDS(
    seurat_obj,
    file.path(
        params$output_dir,
        paste(
            GaitiLabUtils::get_name(params$input_file),
            "reduced_size.rds",
            sep = "_"
        )
    )
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
