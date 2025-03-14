#!/usr/local/bin/_entrypoint.sh Rscript

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)


# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
# Comment/uncomment depending on whether you have an internal package based on the 'R' directory created with usethis::create_package()
# devtools::load_all("./", export_all = FALSE)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Determine consensus for the CCI results of a sample",
    default_output_dir = "output/400_consensus_and_RRA",
    default_log_file = NULL,
    default_log_dir = "output"
)
parser$add_argument(
    "-a",
    "--alpha",
    type = "numeric",
    default = 0.05,
    help = "Significance threshold"
)
parser$add_argument(
    "-id",
    "--sample_id",
    type = "character",
    default = 1,
    help = "Sample id"
)
parser$add_argument(
    "-n",
    "--n_perm",
    type = "numeric",
    default = 1000,
    help = "Number of permutations"
)
parser$add_argument(
    "--cellchat_obj",
    type = "character",
    default = "",
    help = "CellChat object"
)
parser$add_argument(
    "--liana_obj",
    type = "character",
    default = "",
    help = "LIANA object"
)
parser$add_argument(
    "--cell2cell_obj",
    type = "character",
    default = "",
    help = "Cell2Cell object"
)
parser$add_argument(
    "--cpdb_obj",
    type = "character",
    default = "",
    help = "CPDB object"
)
parser$add_argument(
    "--run_dir",
    type = "character",
    default = "",
    help = "Path to run directory"
)
parser$add_argument(
    "--task_id",
    type = "character",
    help = "Task ID from Nextflow, only needed in Nextflow pipeline",
    default = NULL
)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
}

create_dir(params$output_dir)
# Set up logging
logr <- init_logging(log_level = params$log_level, log_file = NULL)
obj_logger <- init_obj_logging(log_file = NULL)
log_info(ifelse(
    interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Parameters:")
log_object(params_ls_to_df(params))


if (
    dir.exists(params$run_dir) &
        ((params$cellchat_obj == "") | (params$liana_obj == "")) |
        (params$cell2cell_obj == "") |
        params$cpdb_obj == ""
) {
    if (params$cellchat_obj == "") {
        params$cellchat_obj <- glue(
            "{params$run_dir}/300_postproc_cellchat/cellchat__{params$sample_id}__postproc.rds"
        )
    }
    if (params$liana_obj == "") {
        params$liana_obj <- glue(
            "{params$run_dir}/301_postproc_liana/liana__{params$sample_id}__postproc.rds"
        )
    }
    if (params$cell2cell_obj == "") {
        params$cell2cell_obj <- glue(
            "{params$run_dir}/302_postproc_cell2cell/cell2cell__{params$sample_id}__postproc.rds"
        )
    }
    if (params$cpdb_obj == "") {
        params$cpdb_obj <- glue(
            "{params$run_dir}/303_postproc_cpdb/cpdb__{params$sample_id}__postproc.rds"
        )
    }
}

log_info("Rank interactions...")
scrnaseq.cellcomm::rra_interactions(
    cellchat_obj = params$cellchat_obj,
    liana_obj = params$liana_obj,
    cell2cell_obj = params$cell2cell_obj,
    cpdb_obj = params$cpdb_obj,
    output_dir = params$output_dir,
    sample_id = params$sample_id,
    n_perm = params$n_perm
)

log_info("Take consensus...")
scrnaseq.cellcomm::take_consensus(
    cellchat_obj = params$cellchat_obj,
    liana_obj = params$liana_obj,
    cell2cell_obj = params$cell2cell_obj,
    cpdb_obj = params$cpdb_obj,
    output_dir = params$output_dir,
    sample_id = params$sample_id,
    alpha = params$alpha
)

log_info("Finished")

log_info("Session Info")
log_object(sessionInfo())

if (!is.null(params$task_id)) {
    write_versions_yml(
        c("scrnaseq.cellcomm", pacman::p_loaded()),
        task_id = params$task_id,
        outdir = params$output_dir
    )
}
