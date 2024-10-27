# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)

# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Combine samples",
        default_output = "output/401_combine_samples"
    )
    parser$add_argument("--input_dir", required = TRUE, help = "Path to input directory")
    parser$add_argument("--metadata", type = "character", help = "Path to metadata (RDS file)", default = NULL)
    parser$add_argument("--sample_var", type = "character", help = "Name of sample variable (default = 'Sample')", default = "Sample")
    parser$add_argument("--condition_var", type = "character", help = "Name of condition variable", default = "Condition_dummy")
    parser$add_argument("--patient_var", type = "character", help = "Name of patient variable", default = "")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- paste0(args$output_dir)
create_dir(output_dir)

log_info("Hard combine files...")
scrnaseq.cellcomm::combine_samples(
    input_dir = args$input_dir,
    metadata = args$metadata,
    patient_var = args$patient_var,
    condition_var = args$condition_var,
    sample_var = args$sample_var,
    output_dir = args$output_dir
)
log_info("Finished!")
