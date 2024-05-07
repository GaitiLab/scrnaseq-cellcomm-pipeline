# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
require(GBMutils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Score pseudobulking",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/pseudobulk_scoring/malignant_only")
    args$input_file <- "output/gbm_regional_study_pseudobulk_malignant_only_by_sample.rds"
    # args$signatures <- "000_misc_local/CIBERSORT_cell_type_signatures.rds"
    # args$signatures <- "000_misc_local/gene_lists/neftel_signatures.rds"
    args$signatures <- "000_misc_local/gene_lists/verhaak_gene_signatures.rds"
    args$suffix <- "verhaak"
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
pacman::p_load(GSVA, readxl, boot, scalop, Hmisc, corrplot)

log_info("Load data...")
gene_exp_mat <- readRDS(args$input_file)

log_info("Load cell type gene signatures...")
# Loading our signatures
signatures_list <- readRDS(args$signatures)

# # https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#3_Overview_of_the_GSVA_functionality
log_info("Compute ssGSEA scores for each signature...")
ssgsea <- GSVA::gsva(as(data.matrix(gene_exp_mat), "dgCMatrix"), signatures_list, method = "ssgsea")

log_info("Compute scalop scores for each signature...")
scalop_scores <- sigScores(data.matrix(gene_exp_mat), sigs = signatures_list)

log_info("Save...")
saveRDS(ssgsea, file = glue("{args$output_dir}/pseudobulking_scored_ssgsea_{args$suffix}.rds"))
saveRDS(scalop_scores, file = glue("{args$output_dir}/pseudobulking_scored_scalop_{args$suffix}.rds"))

log_info("COMPLETED!")
