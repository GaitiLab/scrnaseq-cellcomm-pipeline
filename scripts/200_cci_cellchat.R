if (dir.exists("/opt/.renv")) {
    print("running with Docker/Singularity")
    renv::load("/opt/.renv")
}
print(renv::paths$library())

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
        description = "Inferring CCIs using CellChat"
    )
    parser$add_argument("-r", "--resource",
        type = "character", help = "Path to custom database with interactions."
    )
    parser$add_argument("-i", "--ident_col",
        type = "character", help = "Column name for cell identity",
        default = "cell_type"
    )
    parser$add_argument("-n", "--n_perm",
        type = "integer", help = "Number of permutations",
        default = 1000L
    )
    parser$add_argument("-id", "--id",
        type = "integer",
        default = NULL, help = "Only required when doing downsampling"
    )
    parser$add_argument("-g", "--gene_expr",
        type = "character",
        default = NULL, help = "Name of RDS file"
    )
    parser$add_argument("-s", "--sampling_grid",
        type = "character",
        default = NULL, help = "Name of RDS file, only required when doing downsampling"
    )
    parser$add_argument("--n_cores",
        default = 1, type = "integer",
        help = "Number of cores to use for parallelization"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/test")
    args$ident_col <- "CellClass_L4"
    args$n_perm <- 100
    args$resource <- glue("{here::here()}/001_data/interactions_db_v2/cellchat_db.rds")
    args$gene_expr <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/100_preprocessing/seurat/6234_2895153_B.rds")
    args$sampling_grid <- NULL
    args$id <- NULL
    args$n_cores <- 1
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
pacman::p_load(reshape2, pbapply)
pacman::p_load_gh("jinworks/CellChat")
options(stringsAsFactors = FALSE)

log_info("Load data...")
seurat_obj <- readRDS(args$gene_expr)

# ---- Run CellChat ----
if (is.null(args$sampling_grid)) {
    log_info("Single run...")
    infer_cellchat(seurat_obj = seurat_obj, output_dir = output_dir, args = args)
} else {
    log_info("Load sampling grid...")
    sampling_grid <- readRDS(args$sampling_grid)

    log_info("Sampling cells...")
    # Subset to sampling grid
    cell_ids <- sampling_grid %>%
        filter(run == args$id) %>%
        select(-run) %>%
        as.matrix() %>%
        as.vector()
    seurat_obj_subset <- subset(seurat_obj, cells = cell_ids)
    log_info("Infer interactions...")
    infer_cellchat(i = args$id, seurat_obj = seurat_obj_subset, output_dir = output_dir, args = args)
}
log_info("COMPLETED!")

devtools::session_info()
