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
        description = "Inferring CCIs using LIANA",
    )
    parser$add_argument("-p", "--n_perm",
        type = "integer",
        default = 1000,
        help = "Number of permutations for permutation testing"
    )
    parser$add_argument("-db", "--interactions_db",
        type = "character", default = "",
        help = "Path to custom database with interactions (RDS)"
    )
    parser$add_argument("-a", "--annot",
        type = "character", default = "cell_type",
        help = "Column in metadata containing the cell type labels"
    )
    parser$add_argument("-g", "--gene_expr",
        type = "character",
        default = NULL, help = "Seurat object with the gene expression (RDS file)"
    )
    parser$add_argument("-n", "--min_cells",
        type = "integer", default = 5, help = "Minimum number of cells required in each cell group for cell-cell communication"
    )
    parser$add_argument("-mp", "--min_pct",
        type = "numeric",
        default = 0.1, help = "Minimum percentage of cells expressing a gene"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/test_dowsampling_implementation/201_cci_liana")
    args$annot <- "CCI_CellClass_L1"
    args$n_perm <- 10
    args$interactions_db <- glue("{here::here()}/data/interactions_db/liana_db.rds")
    args$gene_expr <- glue("{here::here()}/output/test_pipeline/100_preprocessing/seurat/6419_cortex__run__1.rds")
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

# Load additional libraries
log_info("Loading libraries...")
pacman::p_load(OmnipathR)
pacman::p_load_gh("saezlab/liana")

# ---- Constants ----
# TODO add "cytotalk" later after finishing nextflow pipeline
methods <- c("natmi", "connectome", "logfc", "sca", "cytotalk")
supp_columns <- c("ligand.expr", "receptor.expr")
permutation_params <- list(
    nperms = args$n_perm
)
assay <- "RNA"

# ---- Loading data ----
log_info("Loading Seurat object...")
seurat_obj <- readRDS(args$gene_expr)

log_info("Loading database with interactions...")
custom_resource <- readRDS(args$interactions_db)

# ---- Run LIANA ----
liana_obj <- liana_wrap(seurat_obj,
    method = methods,
    resource = "custom",
    external_resource = custom_resource,
    idents_col = args$annot,
    supp_columns = supp_columns,
    return_all = TRUE,
    permutation.params = permutation_params,
    assay = assay,
    min_cells = args$min_cells,
    expr_prop = args$min_pct
)
log_info("Save LIANA results...")
saveRDS(liana_obj,
    file = glue("{args$output_dir}/liana__{get_name(args$gene_expr)}.rds")
)

log_info("COMPLETED!")
