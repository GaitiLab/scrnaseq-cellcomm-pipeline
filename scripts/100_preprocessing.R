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
    parser$add_argument("--min_cells",
        type = "numeric",
        default = 5, help = "Minimum number of cells per annotation"
    )
    parser$add_argument("--interactions_db",
        type = "character",
        default = NULL, help = "Path to interactions db"
    )
    parser$add_argument("--celltypes_oi",
        type = "character",
        default = NULL, help = "Path to cell types of interest"
    )
    parser$add_argument("-fn", "--first_n", type = "numeric", default = 0, help = "First N cell types that have to be present from given celltypes_oi file")
    parser$add_argument("--min_cell_types", type = "numeric", default = 2, help = "Minimum number of cell types to be present after filtering")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/CellClass_L1/100_preprocessing")
    args$interactions_db <- glue("{here::here()}/data/interactions_db/interactions_ref.rds")
    args$input_file <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/001_data/GBM/split_by_Sample/6237_2222190_F.rds"
    args$annot <- "CellClass_L1"
    args$min_cells <- 200
    # args$celltypes_oi <- glue("{here::here()}/data/celltypes_oi.txt")
    args$celltypes_oi <- ""
    args$first_n <- 0
    args$min_cell_types <- 3
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

if ((!is.null(args$celltypes_oi))) {
    log_info("Check if cell types of interest file exists...")
    if (file.exists(args$celltypes_oi)) {
        log_info("Loading cell types of interest...")
        celltypes_oi <- read.table(args$celltypes_oi, sep = "\t") %>% pull(V1)
    }
} else {
    log_info("No cell types of interest provided, use all celltypes in object...")
    celltypes_oi <- seurat_obj[[args$annot]] %>%
        pull() %>%
        unique()
    log_info("Celltypes of interest: {paste0(celltypes_oi, collapse = ', ')}")
}

if (args$min_cells < min_cells) {
    stop(glue("Minimum number of cells per annotation must be >= {min_cells}"))
}

# Get list of genes from interactions db for filtering
log_info("Loading interactions db...")
interactions_db <- readRDS(args$interactions_db)

log_info("Extract all unique genes from the interactions db...")
genes_oi <- interactions_db %>%
    select(source_genesymbol, target_genesymbol) %>%
    unlist() %>%
    unique() %>%
    str_split(., "_") %>%
    unlist() %>%
    unique()

# genes_oi <- c(str_split(interactions_db$interaction, "__", simplify = TRUE))
# genes_oi <- c(str_split(genes_oi, ":", simplify = TRUE))
# genes_oi <- unique(genes_oi[genes_oi != ""])

common_cell_types <- intersect(unique(seurat_obj@meta.data[[args$annot]]), celltypes_oi)

if (args$first_n == 0) {
    log_info("Use all available cell types...")
    check_first_n <- TRUE
} else {
    log_info(glue("First N={first_n} cell types in given cell types of interest need to be present..."))
    check_first_n <- celltypes_oi[1:first_n] %in% common_cell_types
}

if (check_first_n && length(common_cell_types) >= args$min_cell_types) {
    log_info("Only keep cell types of interest...")
    cells_to_keep <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[args$annot]] %in% celltypes_oi]
    seurat_obj <- subset(seurat_obj,
        cells = cells_to_keep
    )
    log_info("Filtering...")
    if (args$min_cells >= min_cells) {
        seurat_obj <- filtering(
            seurat_obj,
            annot = args$annot,
            min_cells = args$min_cells, genes_oi
        )
    } else {
        log_info("Skipping filtering...")
    }
    log_info(glue(
        "Check number of cell types after filtering (>= {args$min_cell_types})..."
    ))
    n_cell_types <- length(unique(seurat_obj@meta.data[[args$annot]]))
    output_name <- str_split(get_name(args$input_file), "__", simplify = TRUE)[1]
    if (n_cell_types >= args$min_cell_types) {
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
    } else {
        log_info("Not enough cell types of interest...")
    }
} else {
    log_info("Not enough cell types of interest...")
}

log_info("COMPLETED!")
