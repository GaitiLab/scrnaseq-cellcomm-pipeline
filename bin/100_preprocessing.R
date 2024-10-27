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
        description = "Preprocessing for CCIs", default_output = "output/100_preprocessing"
    )
    parser$add_argument("--input_file",
        type = "character",
        default = NULL, help = "Path to Seurat object"
    )
    parser$add_argument("--annot",
        type = "character",
        default = "CellClass_L1", help = "Annotation to use for filtering"
    )
    parser$add_argument("-n", "--min_cells",
        type = "integer", default = 5, help = "Minimum number of cells required in each cell group for cell-cell communication"
    )
    parser$add_argument("--is_confident", type = "numeric", default = 0, help = "Filter confident cells (1) or not (0); only relevant for internal project")
    parser$add_argument("--sample_id", type = "character", default = NULL, help = "Sample ID")

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
options(Seurat.object.assay.version = "v4")

log_info("Loading Seurat object...")
create_dir(args$output_dir)
log_info("Prepare data...")
scrnaseq.cellcomm::prepare_data(
    input_file = args$input_file,
    annot = args$annot,
    output_dir = args$output_dir,
    is_confident = args$is_confident,
    min_cells = args$min_cells,
    sample_id = args$sample_id
)
log_info("Finished!")
