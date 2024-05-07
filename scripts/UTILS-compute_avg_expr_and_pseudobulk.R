# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)

# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)
if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Compute (1) average expression, (2) pseudobulking using Seurat",
    )
    parser$add_argument("--input_file", help = "Seurat input file", type = "character")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
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

log_info("Load seurat object...")
obj <- readRDS(args$input_file)

log_info("Compute Average Expression...")
# GBM - For scatterplot
# Only keep `Confident_Annotation` cells
avg_expr <- AverageExpression(subset(obj, subset = Confident_Annotation),
    assays = "SCT",
    group.by = c("Sample", "CCI_CellClass_L2_2"),
    slot = "counts"
)$SCT

log_info("Save results...")
saveRDS(avg_expr, file = glue("{args$output_dir}/average_expression_by_Sample_CCI_CellClass_L2_2.rds"))

log_info("Pseudobulk by sample...")
# GBM - For Verhaak subtyping
# Including all cells
pseudobulk <- AggregateExpression(object = obj, group.by = "Sample")$RNA

log_info("Save results...")
saveRDS(pseudobulk, file = glue("{args$output_dir}/gbm_regional_study_pseudobulk_by_Sample.rds"))

log_info("COMPLETED!")
