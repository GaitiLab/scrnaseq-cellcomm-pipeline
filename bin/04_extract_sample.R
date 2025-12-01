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
    description = "Extract sample from Seurat object",
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
    help = "Path to input directory"
)
parser$add_argument(
    "--sample_var",
    type = "character",
    default = "sample_id",
    help = "Name of sample variable, necessary for splitting (default='sample_id')"
)
parser$add_argument(
    "--sample_id",
    type = "character",
    help = "Sample ID for subsetting"
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
            "reduced_size_subset.rds",
            sep = "_"
        )
    )
    params$sample_var <- "Sample"
    params$sample_id <- "Sample_2"
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

seurat_obj <- readRDS(params$input_file)
log_info("Loaded Seurat object.")

cell_ids <- seurat_obj@meta.data |>
    dplyr::filter(!!sym(params$sample_var) == params$sample_id) |>
    row.names()
log_info("Extracted cell IDs for sample of interest.")

seurat_obj <- subset(seurat_obj, cells = cell_ids)
log_info("Subsetted Seurat object.")

# TODO move to GaitiLabUtils
get_assay_version <- function(seurat_obj, assay = "RNA") {
    if (
        stringr::str_detect(
            as.character(class(seurat_obj[[assay]])[1]),
            as.character(5)
        )
    ) {
        return("v5")
    }
    return("v3/v4")
}
assay_version <- seurat_obj |> get_assay_version()
log_info(paste("Detected Seurat assay version:", assay_version))

if (assay_version == "v5") {
    # Recreate as normalization requires a v4 assay (see options() at start of script)
    seurat_obj_recreated <- Seurat::CreateSeuratObject(
        counts = Seurat::GetAssayData(
            seurat_obj,
            layer = "counts",
            assay = "RNA"
        ),
        meta.data = seurat_obj@meta.data
    )
    log_info("Recreated Seurat object with assay version v4")

    saveRDS(
        seurat_obj_recreated,
        file.path(params$output_dir, paste0(params$sample_id, ".rds"))
    )
    log_info("Saved recreated Seurat object.")
} else {
    saveRDS(
        seurat_obj,
        file.path(params$output_dir, paste0(params$sample_id, ".rds"))
    )
    log_info("Saved Seurat object.")
}

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
