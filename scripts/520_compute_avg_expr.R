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
        description = "Get metadata",
    )
    parser$add_argument("--gene_exp", type = "character", help = "Gene expression file")
    parser$add_argument("--interactions_db", type = "character", help = "Interactions database")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
    args$interactions_db <- "001_data_local/interactions_db/interactions_ref.rds"
    args$gene_exp <- glue("{here::here()}/output/CellClass_L4_min3_types/100_preprocessing/seurat/6419_enhancing_border.rds")
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
pacman::p_load(Seurat, ggplot2, ggrepel, ggpubr)

log_info("Load interactions...")
interactions_db <- readRDS(args$interactions_db)
ligands <- interactions_db %>%
    pull(genename_a) %>%
    unique()
receptors <- interactions_db %>%
    pull(genename_b) %>%
    unique()

gene_exp <- readRDS(args$gene_exp)

genes_oi <- c(ligands, receptors)

log_info("Determine average expression per cell type...")
avg_expr <- Seurat::DotPlot(gene_exp, features = intersect(rownames(gene_exp), genes_oi), group.by = "CellClass_L4")$data
rownames(avg_expr) <- NULL

log_info("Writing output...")
saveRDS(avg_expr, file = glue("{args$output_dir}/{get_name(args$gene_exp)
}.rds"))

log_info("COMPLETED!")
