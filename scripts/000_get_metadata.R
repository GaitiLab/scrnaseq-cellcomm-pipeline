# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# require(GBMutils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Get metadata from Seurat object",
    )
    parser$add_argument("-i", "--input_file",
        type = "character",
        default = NULL, help = "Path to Seurat object"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$input_file <- glue("{here::here()}/data/test_seurat_obj_10x.rds")
    args$output_dir <- glue("{here::here()}/test_output/000_misc")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- paste0(args$output_dir)
create_dir(output_dir)

# Load additional libraries
pacman::p_load(Seurat)

log_info("Load seurat_object")
seurat_obj <- readRDS(args$input_file)

# REMARK: Only relevant for GBM project.
# log_info(glue("Number of cells in object, BEFORE: {ncol(seurat_obj)}"))
# seurat_obj <- subset(seurat_obj, subset = Confident_Annotation == TRUE)
# log_info(glue("Number of cells in object, After: {ncol(seurat_obj)}"))

log_info("Saving metadata...")
saveRDS(seurat_obj@meta.data, glue("{output_dir}/{get_name(args$input_file)}__metadata.rds"))
write.csv(seurat_obj@meta.data, glue("{output_dir}/{get_name(args$input_file)}__metadata.csv"))

log_info("COMPLETED!")
