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
        description = "Split integrated objects into individual samples", default_output = "output/000_data/split_by_Sample"
    )
    parser$add_argument("-i", "--input_file",
        type = "character",
        default = NULL, help = "Path to input directory"
    )
    parser$add_argument("--sample_var", type = "character", default = "Sample", help = "Name of sample variable, necessary for splitting")
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

# Load additional libraries
options(Seurat.object.assay.version = "v4")

log_info("Split Seurat object into samples...")
scrnaseq.cellcomm::split_seurat_into_samples(
    input_file = args$input_file,
    output_dir = args$output_dir,
    sample_var = args$sample_var
)
log_info("Finished!")
