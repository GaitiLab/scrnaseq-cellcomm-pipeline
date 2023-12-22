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
        description = "Create run grid",
    )
    parser$add_argument("-i", "--input_dir",
        type = "character", default = NULL,
        help = "Path to directory containing Seurat objects"
    )
    parser$add_argument("-n", "--nruns",
        type = "integer", default = 500,
        help = "Number of runs to perform"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/parsebio_6245/000_misc")
    args$nruns <- 500
    args$input_dir <- glue("{here::here()}/output/parsebio_6245/100_preprocessing/seurat")
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
samples <- list.files(args$input_dir)

# Create run grid
grid <- data.frame(
    list(
        id = rep(seq_len(args$nruns), times = length(samples)),
        sample_id = rep(samples, each = args$nruns)
    )
)
write.csv(grid, glue("{args$output_dir}/grid.csv"), row.names = FALSE)
log_info(glue("Number of runs: {nrow(grid)}"))
log_info("COMPLETED!")
