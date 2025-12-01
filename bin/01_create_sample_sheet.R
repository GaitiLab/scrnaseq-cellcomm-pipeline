#!/usr/local/bin/_entrypoint.sh Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
options(Seurat.object.assay.version = "v4")
pacman::p_load(GaitiLabUtils, glue, data.table, tidyverse, stringr, Seurat)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Create sample sheet",
    default_output_dir = file.path(
        "output",
        "cci_pipeline",
        "01_prepare_data",
        "01_misc"
    )
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
    default = "sample_id",
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

params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    params$input_file <- file.path(
        "output",
        "cci_pipeline",
        "01_prepare_data",
        "01_misc",
        paste("example_data", "metadata.rds", sep = "__")
    )
    params$annot <- "seurat_annotations"
    params$min_cells <- 70
    params$sample_var <- "Sample"
}

create_dir(params$output_dir)
# Set up logging
log_file <- set_logfile(params)
logr <- init_logging(log_level = params$log_level, log_file = log_file)
obj_logger <- init_obj_logging(log_file = log_file)
log_info(ifelse(
    interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Parameters:")
log_object(params_ls_to_df(params))

# ---- Check arguments ----
# checked_path <- is_valid_path(
#     params$input_file,
#     required_file_extension = "rds"
# )
# if (!checked_path) {
#     stop("Given input file is not a valid path.")
# }

# ---- Workflow ----
meta_df <- readRDS(params$input_file)
log_info("Loaded dataframe.")

n_cells_by_sample_id_and_celltype_df_long <- meta_df |>
    # Duckdb cannot handle 'factor' types
    dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.character)) |>
    duckplyr::as_duckdb_tibble() |>
    dplyr::group_by(
        !!dplyr::sym(params$sample_var),
        !!dplyr::sym(params$annot)
    ) |>
    dplyr::count(name = "n_cells") |>
    ungroup()
log_info("Computed no. cells for each sample ID x cell type pair.")

# Ensure that only samples that have enough cells for at least 2 cell types to infer interactions
passing_samples <- n_cells_by_sample_id_and_celltype_df_long |>
    dplyr::mutate(has_enough_cells = n_cells >= params$min_cells) |>
    dplyr::group_by(!!dplyr::sym(params$sample_var)) |>
    dplyr::summarise(n_celltypes_with_enough_cells = sum(has_enough_cells)) |>
    dplyr::ungroup() |>
    # Only keep samples that have at least 2 cell types with enough cells
    dplyr::filter(n_celltypes_with_enough_cells >= 2) |>
    dplyr::pull(!!dplyr::sym(params$sample_var))
log_info("Extracted sample IDs.")

# Dataframe with no. cells for each cell type
samplesheet <- n_cells_by_sample_id_and_celltype_df_long |>
    dplyr::filter(!!dplyr::sym(params$sample_var) %in% passing_samples) |>
    tidyr::pivot_wider(
        names_from = !!dplyr::sym(params$annot),
        values_from = n_cells,
        values_fill = 0
    ) |>
    dplyr::rename(sample_id = !!sym(params$sample_var))
log_info("Created samplesheet.")

write.csv(
    samplesheet,
    file = file.path(params$output_dir, "samplesheet.csv"),
    quote = FALSE,
    row.names = FALSE
)
log_info("Save samplesheet as csv.")

log_info("Finished")

log_info("Session Info")
log_object(sessionInfo())

if (!is.null(params$nf_process_id)) {
    write_versions_yml(
        c("scrnaseq.cellcomm", pacman::p_loaded()),
        task_id = params$nf_process_id,
        outdir = params$output_dir
    )
}
