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
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/000_misc_local")
    args$sample_sheet <- glue("{here::here()}/000_misc_local/gbm_regional_study__20231115.xlsx")
    args$sheet_name <- "CellClass_L4"
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
pacman::p_load(readxl)

log_info("Read sample sheet...")
excel_file <- read_excel(args$sample_sheet, sheet = args$sheet_name)

log_info("Filter...")
excel_file_filtered <- excel_file %>% filter(`N cell types` >= 2, Region != "NC")

log_info("Formatting...")
out <- excel_file_filtered %>%
    rownames_to_column() %>%
    select(rowname, Sample)

log_info("Save output...")
write.table(out, file = glue("{args$output_dir}/samples_oi.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)

log_info("COMPLETED!")
