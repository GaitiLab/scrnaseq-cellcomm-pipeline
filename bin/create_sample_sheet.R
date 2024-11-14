#!/usr/bin/env Rscript

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Create sample sheet", default_output = "."
    )
    parser$add_argument("-i", "--input_file",
        type = "character",
        default = NULL, help = "Path to input directory"
    )
    parser$add_argument("--sample_var", type = "character", default = "Sample", help = "Name of sample variable, necessary for splitting")
    parser$add_argument("--annot",
        type = "character",
        default = "CellClass_L1", help = "Annotation to use for filtering"
    )
    parser$add_argument("-n", "--min_cells",
        type = "integer", default = 5, help = "Minimum number of cells required in each cell group for cell-cell communication"
    )
    parser$add_argument("--is_confident", type = "numeric", default = 0, help = "Filter confident cells (1) or not (0); only relevant for internal project")

    params <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    params <- list()
    params$log_level <- 5
    params$output_dir <- glue("{here::here()}/output/")
    params$input_file <- "output/000_data/example_data__metadata.rds"
    params$annot <- "seurat_annotations"
    params$sample_var <- "Sample"
    params$min_cells <- 75
    params$is_confident <- 0
}

# Set up logging
logr <- init_logging(log_level = params$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(params$output_dir)

# Load additional libraries

log_info("Load data...")
meta_df <- readRDS(params$input_file)

if (params$is_confident) {
    log_info("Filter on confident annotation ('Confident_Annotation')...")
    # NOTE only used for internal GBM project
    meta_df <- meta_df %>% filter(Confident_Annotation)
}

n_by_celltype <- meta_df %>%
    group_by(!!sym(params$sample_var), !!sym(params$annot)) %>%
    count() %>%
    pivot_wider(names_from = !!sym(params$annot), values_from = n, values_fill = 0)

passing_samples <- n_by_celltype %>%
    summarise(across(where(is.numeric), ~ .x > params$min_cells)) %>%
    rowwise(!!sym(params$sample_var)) %>%
    # For each sample, Count number of cell types with n > min_cells
    summarise(n_cell_types = sum(c_across(where(is.logical)))) %>%
    # Only keep samples that have at least 2 cell types with enough cells
    filter(
        n_cell_types >= 2
    ) %>%
    pull(!!sym(params$sample_var))


log_info("Extract and save final sample sheet...")
write.csv(
    n_by_celltype %>%
        filter(!!sym(params$sample_var) %in% passing_samples) %>%
        rename(Sample = !!sym(params$sample_var)),
    file = file.path(params$output_dir, "sample_sheet.csv"), quote = FALSE, row.names = FALSE
)

log_info("Finished!")
