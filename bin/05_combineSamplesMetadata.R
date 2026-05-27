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
    "--input_dir",
    type = "character",
    default = NULL,
    help = "Path to input directory"
)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    params$input_dir <- file.path(
        here::here(),
        "internal",
        "stromaProject",
        "byDonorId",
        "output",
        "01_prepare_data",
        "02_intermediate_objects"
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

.LoadData <- function(path) {
    log_info("Path: ", path)
    return(readRDS(path)@meta.data)
}

# ---- Workflow ----
seuratObjPaths <- list.files(
    params$input_dir,
    pattern = ".rds",
    full.names = TRUE
)
log_info("No. Seurat objects: ", length(seuratObjPaths))

metaDf <- seuratObjPaths |> map_dfr(.LoadData)
log_info("Combined metadata from all samples.")


saveRDS(metaDf, file.path(params$output_dir, "metadata.rds"))
log_info("Saved metadata as rds.")

write.csv(
    metaDf,
    file.path(params$output_dir, "metadata.csv")
)
log_info("Saved metadata csv.")

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
