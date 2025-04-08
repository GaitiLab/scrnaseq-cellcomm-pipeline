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
    description = "Extract sample from Seurat object",
    default_output_dir = "output",
    default_log_file = NULL,
    default_log_dir = "output"
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
    default = "Sample",
    help = "Name of sample variable, necessary for splitting (default='Sample')"
)
parser$add_argument(
    "--sample_id",
    type = "character",
    help = "Sample ID for subsetting"
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

# Load additional libraries
options(Seurat.object.assay.version = "v4")

pacman::p_load(Seurat)

log_info("Load Seurat object...")
seurat_obj <- readRDS(params$input_file)

log_info("Extract cell IDs for sample of interest...")
cell_ids <- seurat_obj@meta.data %>%
    filter(!!sym(params$sample_var) == params$sample_id) %>%
    row.names(.)

log_info("Subset Seurat object...")
seurat_obj <- subset(seurat_obj, cells = cell_ids)

log_info("Check Seurat assay version...")
assay_version <- as.character(class(seurat_obj[["RNA"]])[1])
if (str_detect(assay_version, as.character(5))) {
    log_info(glue("Seurat assay version: {assay_version}"))
    log_info("Recreate Seurat object...")
    # Needed for normalization
    seurat_obj_recreated <- CreateSeuratObject(
        counts = GetAssayData(seurat_obj, layer = "counts", assay = "RNA"),
        meta.data = seurat_obj@meta.data
    )
    saveRDS(
        seurat_obj_recreated,
        file.path(params$output_dir, paste0(params$sample_id, ".rds"))
    )
} else {
    saveRDS(
        seurat_obj,
        file.path(params$output_dir, paste0(params$sample_id, ".rds"))
    )
}

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
