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
        description = "Determine consensus for a sample",
        default_output = "output/400_consensus_and_RRA"
    )
    parser$add_argument("-a", "--alpha",
        type = "numeric",
        default = 0.05, help = "Significance threshold"
    )
    parser$add_argument("-id", "--sample_id",
        type = "character", default = 1,
        help = "Sample id"
    )
    parser$add_argument("-n", "--n_perm",
        type = "numeric", default = 1000,
        help = "Number of permutations"
    )
    parser$add_argument("--cellchat_obj",
        type = "character", default = "",
        help = "CellChat object"
    )
    parser$add_argument("--liana_obj",
        type = "character", default = "",
        help = "LIANA object"
    )
    parser$add_argument("--cell2cell_obj",
        type = "character", default = "",
        help = "Cell2Cell object"
    )
    parser$add_argument("--cpdb_obj",
        type = "character", default = "",
        help = "CPDB object"
    )
    parser$add_argument("--run_dir",
        type = "character", default = "",
        help = "Path to run directory"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$run_dir <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/GBM_CCI_Analysis/output/CCI_CellClass_L2_2_reassigned_samples_confident_only_FINAL"
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/test_individual_scripts/400_consensus_and_RRA")
    args$sample_id <- "6234_2895153_A"
    args$alpha <- 0.05
    args$cellchat_obj <- glue("{args$run_dir}/300_postproc_cellchat/cellchat__{args$sample_id}__postproc.rds")
    args$liana_obj <- glue("{args$run_dir}/301_postproc_liana/liana__{args$sample_id}__postproc.rds")
    args$cell2cell_obj <- glue("{args$run_dir}/302_postproc_cell2cell/cell2cell__{args$sample_id}__postproc.rds")
    args$cpdb_obj <- glue("{args$run_dir}/303_postproc_cpdb/cpdb__{args$sample_id}__postproc.rds")
    args$n_perm <- 1000
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

if (dir.exists(args$run_dir) & ((args$cellchat_obj == "") | (args$liana_obj == "")) | (args$cell2cell_obj == "") | args$cpdb_obj == "") {
    if (args$cellchat_obj == "") {
        args$cellchat_obj <- glue("{args$run_dir}/300_postproc_cellchat/cellchat__{args$sample_id}__postproc.rds")
    }
    if (args$liana_obj == "") {
        args$liana_obj <- glue("{args$run_dir}/301_postproc_liana/liana__{args$sample_id}__postproc.rds")
    }
    if (args$cell2cell_obj == "") {
        args$cell2cell_obj <- glue("{args$run_dir}/302_postproc_cell2cell/cell2cell__{args$sample_id}__postproc.rds")
    }
    if (args$cpdb_obj == "") {
        args$cpdb_obj <- glue("{args$run_dir}/303_postproc_cpdb/cpdb__{args$sample_id}__postproc.rds")
    }
}

log_info("Rank interactions...")
scrnaseq.cellcomm::rra_interactions(
    cellchat_obj = args$cellchat_obj,
    liana_obj = args$liana_obj,
    cell2cell_obj = args$cell2cell_obj,
    cpdb_obj = args$cpdb_obj,
    output_dir = args$output_dir,
    sample_id = args$sample_id,
    n_perm = args$n_perm
)

log_info("Take consensus...")
scrnaseq.cellcomm::take_consensus(
    cellchat_obj = args$cellchat_obj,
    liana_obj = args$liana_obj,
    cell2cell_obj = args$cell2cell_obj,
    cpdb_obj = args$cpdb_obj,
    output_dir = args$output_dir,
    sample_id = args$sample_id,
    alpha = args$alpha
)
log_info("Finished!")
