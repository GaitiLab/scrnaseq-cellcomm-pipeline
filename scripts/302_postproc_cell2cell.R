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
        description = "Post-process Cell2Cell",
        default_output = "output/302_postproc_cell2cell"
    )
    parser$add_argument("-is", "--input_interactions_scores",
        type = "character", default = NULL,
        help = "CSV file with interaction scores from cell2cell"
    )
    parser$add_argument("-ip", "--input_interactions_pval",
        type = "character", default = NULL,
        help = "CSV file with pvalues from cell2cell"
    )
    parser$add_argument("-s", "--sample_id",
        type = "character", default = NULL,
        help = "Sample ID"
    )
    parser$add_argument("--ref_db", type = "character", default = NULL, help = "Path to interactions database")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- "output/test_individual_scripts/302_postproc_cell2cell"
    args$input_interactions <- "output/test_pipeline/202_cci_cell2cell/cell2cell__Sample_6.csv"
    args$sample_id <- "Sample_6"
    args$ref_db <- "data/interactions_db/ref_db.rds"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

log_info("Standardize format of Cell2Cell results...")
scrnaseq.cellcomm::format_cell2cell_wrapper(
    input_interactions_scores = args$input_interactions_scores,
    input_interactions_pval = args$input_interactions_pval,
    output_dir = args$output_dir,
    sample_id = args$sample_id,
    ref_db = args$ref_db
)

log_info("Finished!")
