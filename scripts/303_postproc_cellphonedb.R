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
        description = "Post-processing CellPhoneDB output",
        default_output = "output/303_postproc_cpdb"
    )
    parser$add_argument("--sample_id",
        default = "",
        type = "character", help = "Sample ID"
    )
    parser$add_argument("--interaction_scores",
        default = "",
        type = "character", help = "Path to CellPhoneDB interaction scores"
    )
    parser$add_argument("--pval",
        default = "",
        type = "character", help = "Path to CellPhoneDB p-values"
    )
    parser$add_argument("--sign_means",
        default = "",
        type = "character", help = "Path to CellPhoneDB significant means"
    )
    parser$add_argument("--means",
        default = "",
        type = "character", help = "Path to CellPhoneDB means"
    )
    parser$add_argument("--ref_db", type = "character", default = NULL, help = "Path to interactions database")

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$sample_id <- "Sample_6"
    args$output_dir <- "output/test_individual_scripts/303_postproc_cpdb"
    args$interaction_scores <- glue("output/test_pipeline/203_cci_cpdb/statistical_analysis_interaction_scores__{args$sample_id}.txt")
    args$pval <- glue("output/test_pipeline/203_cci_cpdb/statistical_analysis_pvalues__{args$sample_id}.txt")
    args$sign_means <- glue("output/test_pipeline/203_cci_cpdb/statistical_analysis_significant_means__{args$sample_id}.txt")
    args$means <- glue("output/test_pipeline/203_cci_cpdb/statistical_analysis_means__{args$sample_id}.txt")
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

log_info("Standardize format of CPDB results...")
scrnaseq.cellcomm::format_cpdb(
    interaction_scores = args$interaction_scores,
    pval = args$pval,
    sign_means = args$sign_means,
    means = args$means,
    output_dir = args$output_dir,
    sample_id = args$sample_id,
    ref_db = args$ref_db
)
log_info("Finished!")
