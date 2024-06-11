# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
# devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Get metadata",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- "output/LP_IMM_perSample"
    args$input_file <- "output/LP_IMM_perSample/000_data/LP_IMM_3B13mut_annotation_scvi__metadata.rds"
    args$annot <- "CellClass_L3_LP"
    args$min_cells <- 5
    args$condition_var <- "Mutation"
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
pacman::p_load(xlsx)

log_info("Load metadata...")
meta <- readRDS(args$input_file)

log_info("Count number of cells by label by sample...")
n_cells_by_label_by_sample <- meta %>%
    # filter(Confident_Annotation) %>%
    group_by(!!sym(args$condition_var), Sample, !!sym(args$annot)) %>%
    count() %>%
    pivot_wider(names_from = !!sym(args$annot), values_from = n, values_fill = 0) %>%
    data.frame()

n_cells_by_label_by_sample_binarized <- meta %>%
    # filter(Confident_Annotation) %>%
    group_by(!!sym(args$condition_var), Sample, !!sym(args$annot)) %>%
    count() %>%
    mutate(n = n >= args$min_cells) %>%
    pivot_wider(names_from = !!sym(args$annot), values_from = n, values_fill = 0) %>%
    data.frame()

if (file.exists((glue("{args$output_dir}/n_cells_by_sample.xlsx")))) {
    file.remove(glue("{args$output_dir}/n_cells_by_sample.xlsx"))
}
log_info("Save excel file...")
write.xlsx(n_cells_by_label_by_sample, file = glue("{args$output_dir}/n_cells_by_sample.xlsx"), sheetName = "Number of cells")
write.xlsx(n_cells_by_label_by_sample_binarized, file = glue("{args$output_dir}/n_cells_by_sample.xlsx"), sheetName = "thresholded (binarized)", append = TRUE)

log_info("COMPLETED!")
