#!/usr/local/bin/_entrypoint.sh Rscript

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
options(Seurat.object.assay.version = "v4")
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, duckplyr)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Inferring CCIs using LIANA",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "02_run_cci",
        "04_liana",
        "01_raw",
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
    "--interactions_db",
    type = "character",
    default = "data/interactions_db/liana_db.rds",
    help = "Path to custom database with interactions an rds file (default = 'data/interactions_db/liana_db.rds')"
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
    "--gene_expr",
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
    "-mp",
    "--min_pct",
    type = "numeric",
    default = 0.1,
    help = "Minimum percentage of cells expressing a gene (default = 0.1; 10%)"
)

params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
}

create_dir(params$output_dir)
# Set up logging
logr <- init_logging(
    log_level = params$log_level,
    log_file = NULL
)
obj_logger <- init_obj_logging(NULL)
log_info(ifelse(
    interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Parameters:")
log_object(params_ls_to_df(params))

# ---- Check arguments ----
arg_paths <- c(params$gene_expr, params$interactions_db)
checked_filepaths <- data.frame(
    path = arg_paths,
    required_file_extension = rep("rds", 2)
) |>
    purrr::pmap_lgl(GaitiLabUtils::is_valid_path) |>
    setNames(nm = arg_paths)
if (!all(checked_filepaths)) {
    stop(
        "Not all valid paths, please check the following inputs\n",
        paste(names(checked_filepaths)[!checked_filepaths], collapse = "\n")
    )
}

liana_obj <- scrnaseq.cellcomm::RunLIANA(
    gene_expr = params$gene_expr,
    interactions_db = params$interactions_db,
    min_cells = params$min_cells,
    min_pct = params$min_pct,
    n_perm = params$n_perm,
    annot = params$annot
)
log_info("Ran LIANA...")

saveRDS(
    liana_obj,
    file = file.path(
        params$output_dir,
        paste("liana", paste0(params$sample_id, ".rds"), sep = "__")
    )
)
log_info("Saved LIANA object.")

log_info("Finished")

log_info("Session Info")
log_object(sessionInfo())

if (!is.null(params$nf_process_id)) {
    write_versions_yml(
        unique(c("scrnaseq.cellcomm", pacman::p_loaded(), "liana")),
        task_id = params$nf_process_id,
        outdir = params$output_dir
    )
}
