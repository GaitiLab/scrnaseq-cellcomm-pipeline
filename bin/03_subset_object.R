#!/usr/local/bin/_entrypoint.sh Rscript

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
options(Seurat.object.assay.version = "v4")
pacman::p_load(
    GaitiLabUtils,
    glue,
    data.table,
    tidyverse,
    stringr,
    duckplyr,
    Seurat
)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Subset Seurat object",
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

parser$add_argument(
    "--sample_var",
    type = "character",
    default = "sample_id",
    help = "Name of sample variable, necessary for splitting (default='sample_id')"
)

parser$add_argument(
    "--samplesheet",
    type = "character",
    default = NULL,
    help = "Path to sample sheet"
)

params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    params$input_file <- file.path(
        "output",
        "cci_pipeline",
        "01_prepare_data",
        "02_intermediate_objects",
        paste(
            "example_data",
            "reduced_size.rds",
            sep = "_"
        )
    )
    params$samplesheet <- file.path(
        "output",
        "cci_pipeline",
        "01_prepare_data",
        "01_misc",
        "samplesheet.csv"
    )
    params$sample_var <- "Sample"
}

create_dir(params$output_dir)
# Set up logging
log_file <- set_logfile(params)
logr <- init_logging(log_level = params$log_level, log_file = log_file)
obj_logger <- init_obj_logging(log_file = log_file)
log_info(ifelse(
    interactive(),
    "Running interactively.",
    "Running from command line/terminal."
))

log_info("Parameters:")
log_object(params_ls_to_df(params))

# ---- Check arguments ----
# arg_paths <- c(params$input_file, params$samplesheet)
# checked_filepaths <- data.frame(
#     path = arg_paths,
#     required_file_extension = c("rds", "csv")
# ) |>
#     purrr::pmap_lgl(GaitiLabUtils::is_valid_path) |>
#     setNames(nm = arg_paths)
# if (!all(checked_filepaths)) {
#     stop(
#         "Not all valid paths, please check the following inputs\n",
#         paste(names(checked_filepaths)[!checked_filepaths], collapse = "\n")
#     )
# }

# ---- Workflow ----
samplesheet <- read.csv(params$samplesheet)
log_info("Loaded samplesheet.")

seurat_obj <- readRDS(params$input_file)
log_info("Loaded Seurat object.")

cell_ids <- seurat_obj@meta.data |>
    dplyr::filter(
        !!sym(params$sample_var) %in% (samplesheet |> pull(sample_id))
    ) |>
    row.names()
log_info("Extracted cell IDs for samples of interest.")

# To further reduce Seurat object, only keep cells from samples that have sufficient no. cell types with enough cells for each of those cell types
seurat_obj <- subset(seurat_obj, cells = cell_ids)
log_info("Subsetted Seurat object.")

saveRDS(
    seurat_obj,
    file = file.path(
        params$output_dir,
        paste(
            GaitiLabUtils::get_name(params$input_file),
            "subset.rds",
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
