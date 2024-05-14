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
        description = "Preprocessing for CCIs",
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
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/test_downsampling_implementation/100_preprocessing")
    args$input_file <- glue("{here::here()}/output/test_downsampling_implementation/split_by_Sample/6419_cortex__run__3.rds")
    args$annot <- "CellClass_L1"
    args$min_cells <- 5
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Loading Seurat object...")
seurat_obj <- readRDS(args$input_file)

log_info("Create output directory...")
output_seurat <- glue("{args$output_dir}/seurat")
output_mtx <- glue("{args$output_dir}/mtx")
create_dir(output_seurat)
create_dir(output_mtx)

# Load additional libraries
log_info("Loading additional libraries...")
pacman::p_load(Seurat, DropletUtils)

# ---- Constants ----
# Need at least 2 cell types for communication
min_cells <- 5

log_info("Loading Seurat object...")
seurat_obj <- readRDS(args$input_file)

# TODO remove before commit (only for GBM project)
if (args$is_confident) {
    seurat_obj <- subset(seurat_obj, subset = Confident_Annotation)
}
log_info(glue("Only keep cell type groups with at least {args$min_cells} cells"))
if (args$min_cells >= min_cells) {
    seurat_obj <- filtering(
        seurat_obj,
        annot = args$annot,
        min_cells = args$min_cells
    )
}

log_info(glue(
    "Check number of cell types after filtering (>= {args$min_cell_types})..."
))

n_cell_types <- length(unique(seurat_obj@meta.data[[args$annot]]))

output_name <- str_split(
    get_name(args$input_file), "__",
    simplify = TRUE
)

if (n_cell_types > 1) {
    log_info("Normalizing data...")
    seurat_obj <- NormalizeData(seurat_obj)

    log_info("Saving Seurat object...")
    saveRDS(
        seurat_obj,
        glue("{output_seurat}/{output_name}.rds")
    )
    log_info("Convert to mtx format (for Cell2Cell)...")
    mat <- seurat_obj[["RNA"]]@data
    write10xCounts(
        glue("{output_mtx}/{output_name}"),
        mat
    )
}
log_info("COMPLETED!")
