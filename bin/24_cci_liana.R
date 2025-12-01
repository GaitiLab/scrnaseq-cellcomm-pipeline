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
    "--interactions_db",
    type = "character",
    default = file.path("data", "interactions_db", "liana_db.rds"),
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
parser$add_argument("--sample_id", type = "character", default = NULL)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    params$gene_expr <- "internal/stromaProject/bySourceSampleId/output/01_prepare_data/03_preprocessed_objects/seurat/50y.rds"
    params$min_cells <- 100

    params$interactions_db <- "assets/interactions_db/liana_db.rds"
    params$n_perm <- 10
    params$annot <- "CellClass_L2"
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
# arg_paths <- c(params$gene_expr, params$interactions_db)
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

# Run all method
methods <- c("natmi", "connectome", "logfc", "sca", "cytotalk")
supp_columns <- c("ligand.expr", "receptor.expr")
# Define the no. of permutations for permutation testing when applicable
permutation_params <- list(nperms = params$n_perm)
# ---- Perform sanity checks ----
# Minimum of 5 cells enforced/required by LIANA
if (params$min_cells < 5) {
    stop("Min cells has to be >= 5...")
}
# In documentation of LIANA `min_pct` actually represents a fraction/proportion, not a percentage. Therefore value should not be greater than 1.
if (params$min_pct > 1) {
    stop("min_pct > 1...")
}


# ---- Loading data
message("Loading Seurat object...")
seurat_obj <- readRDS(params$gene_expr)


assay <- "RNA"
if (!assay %in% Seurat::Assays(seurat_obj)) {
    stop("`RNA` assay is not present...")
}
message("Loading database with interactions...")
custom_resource <- readRDS(params$interactions_db)

# ---- Run LIANA
lianaObj <- liana::liana_wrap(
    seurat_obj,
    method = methods,
    resource = "custom",
    external_resource = custom_resource,
    idents_col = params$annot,
    supp_columns = supp_columns,
    return_all = TRUE,
    permutation.params = permutation_params,
    assay = assay,
    min_cells = params$min_cells,
    expr_prop = params$min_pct
)

# lianaObj <- scrnaseq.cellcomm::RunLIANA(
#     gene_expr = params$gene_expr,
#     interactions_db = params$interactions_db,
#     min_cells = params$min_cells,
#     min_pct = params$min_pct,
#     n_perm = params$n_perm,
#     annot = params$annot
# )
# log_info("Ran LIANA...")

saveRDS(
    lianaObj,
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
    GaitiLabUtils::write_versions_yml(
        unique(c("scrnaseq.cellcomm", pacman::p_loaded(), "liana")),
        task_id = params$nf_process_id,
        outdir = params$output_dir
    )
}
