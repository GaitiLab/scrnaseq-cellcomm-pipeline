# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
cmd_args <- commandArgs(trailingOnly = FALSE)
has_script_filepath <- startsWith(cmd_args, "--file=")
if (sum(has_script_filepath)) {
    setwd(dirname(unlist(strsplit(cmd_args[has_script_filepath], "=")))[2])
}

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)
if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Generate downsampling runs"
    )
    parser$add_argument("-sp", "--input_file",
        type = "character",
        default = NULL, help = "Path to Seurat object"
    )
    parser$add_argument("-n", "--nruns",
        type = "integer",
        default = "100", help = "Number of runs to perform"
    )
    parser$add_argument("-a", "--annot",
        type = "character", default = NULL,
        help = "Annotation to use for filtering"
    )
    parser$add_argument("-m", "--n_cells", help = "Number of cells to sample", default = NULL)
    parser$add_argument("-c", "--celltypes_oi", help = "File with cell types of interest", default = NULL)
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$input_file <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/cell-cell-interactions/output/parsebio_restrictive/100_preprocessing/seurat/6245_4972288_D
__Tumour_edge__Batch_3.rds"
    args$nruns <- 1000
    args$annot <- "custom_annot"
    args$output_dir <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/cell-cell-interactions/output/test/"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

if (!file.exists(args$input_file)) {
    stop(glue("Input file {args$input_file} does not exist!"))
}

log_info("Create output directory...")
create_dir(args$output_dir)

# Load additional libraries
pacman::p_load(Seurat, pbapply)

log_info("Load Seurat object...")
seurat_obj <- readRDS(args$input_file)

# Get the cell types of interest
if (!is.null(args$celltypes_oi)) {
    celltypes_oi <- read.table(args$celltypes_oi) %>%
        pull(V1)
} else {
    celltypes_oi <- unique(seurat_obj[[args$annot]])
}

# Check whether the user entered a specific number of cells to sample,
# otherwise base number of cells on the limiting cell type.
if (!is.null(args$n_cells)) {
    log_info(glue("Fix number of cells for sampling to {args$n_cells}"))
    n_cells <- args$n_cells
} else {
    log_info("Determine number of cells to sample for each cell type...")
    n_cells <- min(table(seurat_obj[[args$annot]]))
}

log_info(glue("Sample cells for each cell type for {args$nruns} runs..."))
sampling_grid <- pblapply(seq_len(args$nruns), sample_iteration,
    annot_name = args$annot,
    cell_types = celltypes_oi,
    n_cells = n_cells,
    seurat_obj = seurat_obj
)

log_info("Create one dataframe")
df <- do.call(
    rbind,
    lapply(seq_len(length(sampling_grid)), convert_iteration_matrix_to_df,
        sampling_grid = sampling_grid
    )
)

log_info("Saving sampling grid...")
saveRDS(df, glue("{args$output_dir}/{get_name(args$input_file)}__grid.rds"))
log_info("COMPLETED!")
