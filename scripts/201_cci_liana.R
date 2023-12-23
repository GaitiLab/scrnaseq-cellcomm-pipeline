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
        description = "Inferring CCIs using LIANA",
    )
    parser$add_argument("-p", "--n_perms",
        type = "integer",
        default = 1000,
        help = "Number of permutations for significant interactions."
    )
    parser$add_argument("-r", "--resource",
        type = "character", default = "Consensus",
        help = "Path to custom database with interactions."
    )
    parser$add_argument("-i", "--ident_col",
        type = "character", default = "cell_type",
        help = "Identity to use for analysis."
    )
    parser$add_argument("-id", "--id", type = "integer")
    parser$add_argument("-g", "--gene_expr",
        type = "character",
        default = NULL, help = "Name of RDS file"
    )
    parser$add_argument("-n", "--n_cells",
        type = "integer", default = 5, help = "Number of cells to use"
    )
    parser$add_argument("-s", "--sampling_grid",
        type = "character", default = NULL, help = "Name of RDS file"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/output_test/201_liana")
    args$ident_col <- "CellClass_L4"
    args$n_perm <- 1000
    args$resource <- glue("{here::here()}/001_data_local/interactions_db_v2/liana_db.rds")
    args$gene_expr <- glue("{here::here()}/output/CellClass_L4_min3_types/100_preprocessing/seurat/6419_cortex.rds")
    # args$sampling_grid <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/cell-cell-interactions/output/parsebio/101_downsampling/6245_4972288_E__Tumour_edge__Batch_3__grid.rds"
    args$sampling_grid <- NULL
    args$id <- 1
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
# methods <- c("natmi", "connectome", "logfc", "sca", "cellphonedb", "cytotalk")
methods <- c("natmi", "connectome", "logfc", "sca", "cytotalk")
supp_columns <- c("ligand.expr", "receptor.expr")
return_all <- FALSE
permutation_params <- list(nperms = args$n_perms)
assay <- "RNA"

# ---- Loading data ----
log_info("Loading Seurat object...")
seurat_obj <- readRDS(args$gene_expr)

log_info("Loading database with interactions...")
custom_resource <- readRDS(args$resource)

# ---- Run LIANA ----
if (is.null(args$sampling_grid)) {
    liana_obj <- liana_wrap(seurat_obj,
        method = methods, resource = "custom",
        external_resource = custom_resource,
        idents_col = args$ident_col,
        supp_columns = supp_columns, return_all = return_all,
        permutation.params = permutation_params, assay = assay
    )
    log_info("Save LIANA results...")
    saveRDS(liana_obj,
        file = glue("{args$output_dir}/liana__{get_name(args$gene_expr)}.rds")
    )
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
    liana_obj <- liana_wrap(seurat_obj_subset,
        method = methods, resource = "custom",
        external_resource = custom_resource,
        idents_col = args$ident_col,
        supp_columns = supp_columns,
        return_all = return_all,
        permutation.params = permutation_params,
        assay = assay
    )
    log_info("Save LIANA results...")
    saveRDS(liana_obj,
        file = glue("{args$output_dir}/liana__{get_name(args$gene_expr)}__{args$id}.rds")
    )
}

log_info("COMPLETED!")
