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
    args$output_dir <- glue("{here::here()}/001_data_local/omnipathr_data")
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
pacman::p_load(OmnipathR)

pathways <- kegg_pathway_annotations()
saveRDS(pathways, glue("{args$output_dir}/omnipath_pathways.rds"))


protein_complexes <- import_omnipath_complexes()
saveRDS(protein_complexes, glue("{args$output_dir}/omnipath_protein_complexes.rds"))


ppi <-
    import_omnipath_interactions(
        entity_types = "protein"
    )
saveRDS(ppi, glue("{args$output_dir}/omnipath_ppi.rds"))

icn <- import_intercell_network()
saveRDS(icn, glue("{args$output_dir}/omnipath_icn.rds"))
