#!/usr/local/bin/_entrypoint.sh Rscript

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, Seurat)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Get metadata from Seurat object",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "01_prepare_data",
        "01_misc"
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
checked_path <- is_valid_path(
    params$input_file,
    required_file_extension = "rds"
)
if (!checked_path) {
    stop("Given input file is not a valid path.")
}

# ---- Workflow ----
seurat_obj <- readRDS(params$input_file)
log_info("Loaded Seurat object.")

out_filename <- GaitiLabUtils::get_name(params$input_file)

saveRDS(
    seurat_obj@meta.data,
    file.path(
        params$output_dir,
        paste(out_filename, "metadata.rds", sep = "__")
    )
)
log_info("Saved metadata as rds.")

write.csv(
    seurat_obj@meta.data,
    file.path(
        params$output_dir,
        paste(out_filename, "metadata.csv", sep = "__")
    )
)
log_info("Saved metadata as csv.")

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
