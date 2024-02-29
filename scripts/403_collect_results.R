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
        description = "Create ExcelSheet with interactions from post-filtered data",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("/Users/joankant/Library/CloudStorage/OneDrive-SharedLibraries-UHN/Wu, Yiyan - Spatial_GBM/Analysis/CCI")
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

sheet_names <- c("CCI_CellClass_L1_conf_min50", "CCI_CellClass_L2_conf_min50")
for (sheet_name in sheet_names) {
    log_info(glue("Annotation: {sheet_name}..."))
    log_info("Load data...")
    obj <- readRDS(glue("{here::here()}/output/{sheet_name}/402_aggregation/402_patient_interactions_mvoted_w_filters.rds"))

    log_info("Filtering data...")
    obj_filtered <- obj %>%
        filter(lenient_condition) %>%
        select(Region_Grouped, source_target, complex_interaction, lenient_condition_n_patients, lenient_condition_n_samples, lenient_condition_samples, lenient_condition_patients) %>%
        distinct() %>%
        as.data.frame()

    log_info("Write data to Excel...")
    filename <- glue("{args$output_dir}/interactions_conf_malign_min50.xlsx")
    write.xlsx(obj_filtered, filename,
        sheetName = sheet_name,
        col.names = TRUE, row.names = FALSE, append = TRUE
    )
}
