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
        description = "Inferring CCIs using CellChat",
        default_output = "output/200_cci_cellchat"
    )
    parser$add_argument("-p", "--n_perm",
        type = "integer",
        default = 1000,
        help = "Number of permutations for permutation testing (default = 1000)"
    )
    parser$add_argument("-db", "--interactions_db",
        type = "character", default = "data/interactions_db/cellchat_db.rds",
        help = "Path to custom database with interactions (RDS) (default = 'data/interactions_db/cellchat_db.rds')"
    )
    parser$add_argument("-a", "--annot",
        type = "character", default = "cell_type",
        help = "Column in metadata containing the cell type labels"
    )
    parser$add_argument("-g", "--gene_expr",
        type = "character",
        default =
            "", help = "Seurat object with the gene expression (RDS file)"
    )
    parser$add_argument("-n", "--min_cells",
        type = "integer", default = 5, help = "Minimum number of cells required in each cell group for cell-cell communication (default = 5)"
    )
    parser$add_argument("-nc", "--n_cores",
        default = 1, type = "integer",
        help = "Number of cores to use for parallelization (default =1)"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$gene_expr <- "output/test_individual_scripts/100_preprocessing/seurat/Sample_6.rds"
    args$output_dir <- "output/test_individual_scripts/200_cci_cellchat"
    args$annot <- "seurat_annotations"
    args$min_cells <- 50
    args$interactions_db <- "data/interactions_db/cellchat_db.rds"
    args$n_perm <- 5
    args$n_cores <- 2
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

options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v4")

# Set up for parallelization
future::plan("multisession", workers = args$n_cores)

log_info("Run CellChat...")
scrnaseq.cellcomm::run_cellchat(
    gene_expr = args$gene_expr,
    annot = args$annot,
    interactions_db = args$interactions_db,
    output_dir = args$output_dir,
    min_cells = args$min_cells,
    n_perm = args$n_perm
)
log_info("Finished!")
