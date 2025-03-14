#!/usr/local/bin/_entrypoint.sh Rscript

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)


# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
# Comment/uncomment depending on whether you have an internal package based on the 'R' directory created with usethis::create_package()
# devtools::load_all("./", export_all = FALSE)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Inferring CCIs using LIANA",
    default_output_dir = "output/201_cci_liana",
    default_log_file = NULL,
    default_log_dir = "output"
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
parser$add_argument(
    "--task_id",
    type = "character",
    help = "Task ID from Nextflow, only needed in Nextflow pipeline",
    default = NULL
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

options(Seurat.object.assay.version = "v4")

log_info("Run LIANA...")
scrnaseq.cellcomm::run_liana(
    gene_expr = params$gene_expr,
    interactions_db = params$interactions_db,
    output_dir = params$output_dir,
    min_cells = params$min_cells,
    min_pct = params$min_pct,
    n_perm = params$n_perm,
    annot = params$annot
)

log_info("Finished")

log_info("Session Info")
log_object(sessionInfo())

if (!is.null(params$task_id)) {
    write_versions_yml(
        unique(c("scrnaseq.cellcomm", pacman::p_loaded(), "liana")),
        task_id = params$task_id,
        outdir = params$output_dir
    )
}
