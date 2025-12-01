#!/usr/local/bin/_entrypoint.sh Rscript

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v4")

pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, duckplyr)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Inferring CCIs using CellChat",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "02_run_cci",
        "02_cellchat",
        "01_raw"
    )
)
parser$add_argument(
    "-p",
    "--n_perm",
    type = "integer",
    default = 1000,
    help = "Number of permutations for permutation testing (default = 1000)"
)
parser$add_argument(
    "-db",
    "--interactions_db_path",
    type = "character",
    default = "data/interactions_db/cellchat_db.rds",
    help = "Path to custom database with interactions (RDS) (default = 'data/interactions_db/cellchat_db.rds')"
)
parser$add_argument(
    "-a",
    "--annot",
    type = "character",
    default = "cell_type",
    help = "Column in metadata containing the cell type labels"
)
parser$add_argument(
    "-g",
    "--gene_expr_path",
    type = "character",
    default = "",
    help = "Seurat object with the gene expression (RDS file)"
)
parser$add_argument(
    "-n",
    "--min_cells",
    type = "integer",
    default = 5,
    help = "Minimum number of cells required in each cell group for cell-cell communication (default = 5)"
)
parser$add_argument(
    "-nc",
    "--n_cores",
    default = 1,
    type = "integer",
    help = "Number of cores to use for parallelization (default =1)"
)
parser$add_argument(
    "--sample_id",
    type = "character",
    default = "Character string with sample ID"
)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs

    params$cci_dir <- "/cluster/projects/gaitigroup/Users/Jiaoyi/breast_scrnaseq/07_output/CCI/BRCA2_BSO_Neg_CellClass_L2"

    params$gene_expr_path <- "internal/stromaProject/bySourceSampleId/output/01_prepare_data/03_preprocessed_objects/seurat/50y.rds"
    params$annot <- "CellClass_L2"

    params$n_perm <- 10
    params$interactions_db_path <- "assets/interactions_db/cellchat_db.rds"
    params$min_cells <- 100
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

# ---- Check arguments ----
# arg_paths <- c(params$gene_expr_path, params$interactions_db)
# checked_filepaths <- data.frame(
#     path = arg_paths,
#     required_file_extension = rep("rds", 2)
# ) |>
#     purrr::pmap_lgl(GaitiLabUtils::is_valid_path) |>
#     setNames(nm = arg_paths)
# if (!all(checked_filepaths)) {
#     stop(
#         "Not all valid paths, please check the following inputs\n",
#         paste(names(checked_filepaths)[!checked_filepaths], collapse = "\n")
#     )
# }

# Set up for parallelization
future::plan("multisession", workers = params$n_cores)

cellchat_results <- scrnaseq.cellcomm::RunCellChat(
    gene_expr = params$gene_expr_path,
    annot = params$annot,
    interactions_db = params$interactions_db_path,
    min_cells = params$min_cells,
    n_perm = params$n_perm
)
log_info("Ran CellChat...")

saveRDS(
    cellchat_results,
    file.path(
        params$output_dir,
        paste("cellchat", params$sample_id, "raw_obj.rds", sep = "__")
    )
)
log_info("Saved CellChat object...")

log_info("Finished")

log_info("Session Info")
log_object(sessionInfo())

if (!is.null(params$nf_process_id)) {
    write_versions_yml(
        unique(c("scrnaseq.cellcomm", pacman::p_loaded(), "CellChat")),
        task_id = params$nf_process_id,
        outdir = params$output_dir
    )
}
