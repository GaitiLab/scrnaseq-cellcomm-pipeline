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
        description = "Get metadata",
    )
    parser$add_argument("--meta", type = "character", help = "Metadata file")
    parser$add_argument("--sample_varname", type = "character", help = "column name in metadata representing the sample")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/000_misc_local/")
    args$meta <- "000_misc_local/gbmap_core__metadata.rds"
    args$sample_varname <- "sample"
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
log_info("Load metadata...")
samples <- readRDS(args$meta) %>%
    pull(!!sym(args$sample_varname)) %>%
    unique()

log_info("Save list of samples...")
output_filename <- paste0(str_remove(get_name(args$meta), "__metadata"), "__samples.txt")
write.table(samples, file = glue("{args$output_dir}/{output_filename}"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

log_info("COMPLETED!")
