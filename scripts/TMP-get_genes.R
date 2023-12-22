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
        description = "Get genes",
    )
    parser$add_argument("--input_file", type = "character", required = TRUE, help = "Input file")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
    # args$input_file <- "output/CellClass_L4_min3_types/100_preprocessing/seurat/6419_cortex.rds"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

log_info("Reading input file...")
# obj <- readRDS(args$input_file)

# log_info("Getting captured genes...")
# captured_genes <- rownames(obj)
# log_info(glue("Number of captured genes: {length(captured_genes)}"))

# log_info("Saving captured genes...")
# saveRDS(captured_genes, file = glue("{args$output_dir}/captured_genes.rds"))

# log_info("COMPLETED!")


obj <- readRDS(glue("{here::here()}/001_data_local/interactions_db_v2/liana_db.rds"))

tmp <- obj %>%
    select(source_genesymbol, target_genesymbol) %>%
    unlist() %>%
    unique() %>%
    str_split(., "_") %>%
    unlist() %>%
    unique()
