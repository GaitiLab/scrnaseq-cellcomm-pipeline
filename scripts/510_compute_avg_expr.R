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
        description = "Compute average expression",
    )
    parser$add_argument("--gene_exp", type = "character", help = "Gene expression file")
    parser$add_argument("--ref_db", type = "character", help = "Interactions database")
    parser$add_argument("--annot", type = "character", help = "Annotation to use for filtering")
    parser$add_argument("--sample_id", type = "character", help = "Sample ID")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
    args$ref_db <- glue("{here::here()}/001_data_local/interactions_db_v2/ref_db.rds")
    args$gene_exp <- glue("{here::here()}/output/CCI_CellClass_L1/100_preprocessing/seurat/6419_enhancing_border.rds")
    args$annot <- "CCI_CellClass_L1"
    args$sample_id <- "6419_enhancing_border"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Load additional libraries
pacman::p_load(Seurat)

log_info("Load gene expression...")
gene_exp <- readRDS(args$gene_exp)

log_info("Determine average expression per cell type...")
avg_expr <- Seurat::DotPlot(gene_exp, features = rownames(gene_exp), group.by = args$annot)$data
rownames(avg_expr) <- NULL

log_info("Writing output...")
saveRDS(avg_expr, file = glue("{args$output_dir}/{args$sample_id}.rds"))

log_info("COMPLETED!")
