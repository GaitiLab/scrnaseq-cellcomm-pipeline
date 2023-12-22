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
    args$output_dir <- glue("{here::here()}/000_misc_local/")
    args$input_file <- glue("{here::here()}/001_data_local/gbm_regional_study__metadata.rds")
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
pacman::p_load(Seurat, xlsx)

output_file <- glue("{args$output_dir}/{str_remove(get_name(args$input_file), '__metadata')}__{get_current_date()}.xlsx")

log_info("Check if output file already exists...")
if (file.exists(output_file)) {
    log_info("Output file already exists. Removing...")
    file.remove(output_file)
}

log_info("Read metadata...")
meta <- readRDS(args$input_file)

log_info("Add custom annotation...")
meta <- meta %>% mutate(CellClass_L4 = case_when(
    startsWith(CellClass_L2, "Malignant") ~ "Malignant", TRUE ~ CellClass_L2
))

log_info("CellClass_L1...")
cellclass_L1 <- meta %>%
    group_by(Sample, Region, CellClass_L1) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = CellClass_L1, values_from = n, values_fill = 0) %>%
    data.frame()

write.xlsx(cellclass_L1, file = output_file, row.names = FALSE, append = TRUE, sheetName = "CellClass_L1", col.names = TRUE)

log_info("CellClass_L2...")
cellclass_L2 <- meta %>%
    group_by(Sample, Region, CellClass_L2) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = CellClass_L2, values_from = n, values_fill = 0) %>%
    data.frame()

write.xlsx(cellclass_L2, file = output_file, row.names = FALSE, append = TRUE, sheetName = "CellClass_L2", col.names = TRUE)

log_info("CellClass_L3...")
cellclass_L3 <- meta %>%
    group_by(Sample, Region, CellClass_L3) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = CellClass_L3, values_from = n, values_fill = 0) %>%
    data.frame()

write.xlsx(cellclass_L3, file = output_file, row.names = FALSE, append = TRUE, sheetName = "CellClass_L3", col.names = TRUE)

log_info("CellClass_L4...")
cellclass_L4 <- meta %>%
    group_by(Sample, Region, CellClass_L4) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = CellClass_L4, values_from = n, values_fill = 0) %>%
    data.frame()
write.xlsx(cellclass_L4, file = output_file, row.names = FALSE, append = TRUE, sheetName = "CellClass_L4", col.names = TRUE)

log_info("COMPLETED!")
