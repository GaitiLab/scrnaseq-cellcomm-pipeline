# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)

# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Remove redundant assays from Seurat object", default_output = "output/000_data"
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
options(Seurat.object.assay.version = "v4")
log_info("Reduce Seurat object size...")
scrnaseq.cellcomm::reduce_seurat_obj(
    input_file = args$input_file,
    output_dir = args$output_dir
)
log_info("Finished!")
