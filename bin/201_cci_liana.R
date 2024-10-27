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
        description = "Inferring CCIs using LIANA", default_output = "output/201_cci_liana"
    )
    parser$add_argument("-p", "--n_perm",
        type = "integer",
        default = 1000,
        help = "Number of permutations for permutation testing (default = 1000)"
    )
    parser$add_argument("-db", "--interactions_db",
        type = "character", default = "data/interactions_db/liana_db.rds",
        help = "Path to custom database with interactions an rds file (default = 'data/interactions_db/liana_db.rds')"
    )
    parser$add_argument("-a", "--annot",
        type = "character", default = "cell_type",
        help = "Column in metadata containing the cell type labels"
    )
    parser$add_argument("-g", "--gene_expr",
        type = "character",
        default = "", help = "Seurat object with the gene expression (RDS file)"
    )
    parser$add_argument("-n", "--min_cells",
        type = "integer", default = 5, help = "Minimum number of cells required in each cell group for cell-cell communication (default = 5)"
    )
    parser$add_argument("-mp", "--min_pct",
        type = "numeric",
        default = 0.1, help = "Minimum percentage of cells expressing a gene (default = 0.1; 10%)"
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
output_dir <- glue("{args$output_dir}")
create_dir(output_dir)
options(Seurat.object.assay.version = "v4")

# Load additional libraries
log_info("Loading libraries...")

log_info("Run LIANA...")
scrnaseq.cellcomm::run_liana(
    gene_expr = args$gene_expr,
    interactions_db = args$interactions_db,
    output_dir = args$output_dir,
    min_cells = args$min_cells,
    min_pct = args$min_pct,
    n_perm = args$n_perm, annot = args$annot
)

log_info("Finished!")
