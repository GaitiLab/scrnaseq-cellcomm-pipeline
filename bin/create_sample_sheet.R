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
    description = "Create sample sheet",
    default_output_dir = "output",
    default_log_file = NULL,
    default_log_dir = "output"
)
parser$add_argument(
    "-i",
    "--input_file",
    type = "character",
    default = NULL,
    help = "Path to input directory"
)
parser$add_argument(
    "--sample_var",
    type = "character",
    default = "Sample",
    help = "Name of sample variable, necessary for splitting"
)
parser$add_argument(
    "--annot",
    type = "character",
    default = "CellClass_L1",
    help = "Annotation to use for filtering"
)
parser$add_argument(
    "-n",
    "--min_cells",
    type = "integer",
    default = 5,
    help = "Minimum number of cells required in each cell group for cell-cell communication"
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
    params$input_file <- "test_data/example_data.rds"
    params$output_dir <- "test_data"
    params$annot <- "seurat_annotations"
    params$min_cells <- 70
}

create_dir(params$output_dir)
# Set up logging
logr <- init_logging(
    log_level = params$log_level,
    log_file = NULL
)
obj_logger <- init_obj_logging(log_file = NULL)
log_info(ifelse(
    interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Parameters:")
log_object(params_ls_to_df(params))

options(Seurat.object.assay.version = "v4")

log_info("Load data...")
meta_df <- readRDS(params$input_file)

n_by_celltype <- meta_df %>%
    group_by(!!sym(params$sample_var), !!sym(params$annot)) %>%
    count() %>%
    pivot_wider(
        names_from = !!sym(params$annot),
        values_from = n,
        values_fill = 0
    )

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
    file = file.path(params$output_dir, "samplesheet.csv"),
    quote = FALSE,
    row.names = FALSE
)

log_info("Finished")

log_info("Session Info")
log_object(sessionInfo())

if (!is.null(params$task_id)) {
    write_versions_yml(
        c("scrnaseq.cellcomm", pacman::p_loaded()),
        task_id = params$task_id,
        outdir = params$output_dir
    )
}
