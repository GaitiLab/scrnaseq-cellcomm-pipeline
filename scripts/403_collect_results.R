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
        description = "Create ExcelSheet with interactions from post-filtered data",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("/Users/joankant/Library/CloudStorage/OneDrive-UHN/Spatial_GBM/Analysis/WIP/scRNAseq/CCI")
    args$input_dir <- "output/CCI_CellClass_L2_2_reassigned_samples_confident_only"
    args$output_name <- "CCI_CellClass_L2_2_reassigned_samples_confident_only"
    args$condition_varname <- "Region"
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
require(xlsx)
output_filename <- glue("{args$output_dir}/{args$output_name}.xlsx")

if (file.exists(output_filename)) {
    file.remove(output_filename)
}

obj <- readRDS(glue("{args$input_dir}/402_aggregation/402c_aggregation_integration.rds"))
obj_filtered <- obj %>%
    filter(!is.na(lenient_condition), lenient_condition, pval < 0.05) %>%
    select(!!sym(args$condition_varname), source_target, complex_interaction, lenient_condition_n_patients, lenient_condition_n_samples, lenient_condition_samples, lenient_condition_patients, pval, LIANA_score, CellPhoneDB_score, CellChat_score) %>%
    distinct() %>%
    as.data.frame()

log_info("Write data to Excel...")
write.xlsx(obj_filtered, glue("{args$output_dir}/{args$output_name}.xlsx"),
    col.names = TRUE, row.names = FALSE, append = TRUE
)
