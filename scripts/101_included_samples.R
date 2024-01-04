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
    args$output_dir <- glue("{here::here()}/output/")
    args$meta <- glue("{here::here()}/output/CCI_CellClass_L1/000_data/gbm_regional_study__metadata.rds")
    args$min_cells <- 100
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
pacman::p_load(xlsx)
meta <- readRDS(args$meta)

cells_by_sample_filename <- glue("{args$output_dir}/cells_by_sample.xlsx")
if (file.exists(cells_by_sample_filename)) {
    file.remove(cells_by_sample_filename)
}

cells_by_sample_L1 <- meta %>%
    select(CCI_CellClass_L1, Sample, Region_Grouped) %>%
    count(Region_Grouped, Sample, CCI_CellClass_L1) %>%
    pivot_wider(names_from = "CCI_CellClass_L1", values_from = "n", values_fill = 0) %>%
    data.frame()

cells_by_sample_L2 <- meta %>%
    select(CCI_CellClass_L2, Region_Grouped, Sample) %>%
    count(Region_Grouped, Sample, CCI_CellClass_L2) %>%
    pivot_wider(names_from = "CCI_CellClass_L2", values_from = "n", values_fill = 0) %>%
    data.frame()


write.xlsx(cells_by_sample_L1,
    file = cells_by_sample_filename,
    sheetName = "CCI_CellClass_L1", append = FALSE
)

write.xlsx(cells_by_sample_L2,
    file = cells_by_sample_filename,
    sheetName = "CCI_CellClass_L2", append = TRUE
)
