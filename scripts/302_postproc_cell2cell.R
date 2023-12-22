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
        description = "Post-process Cell2Cell",
    )
    parser$add_argument("-i", "--input_interactions",
        type = "character", default = NULL,
        help = "Directory or file with Cell2Cell results"
    )
    parser$add_argument("-s", "--sample_id",
        type = "character", default = NULL,
        help = "Sample ID"
    )
    parser$add_argument("--ref_db",
        type = "character", default = NULL,
        help = "Path to reference database"
    )

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types/302_postproc_cell2cell")
    args$input_interactions <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/202_cci_cell2cell/cell2cell__6234_2895153_A.csv")
    args$sample_id <- "6234_2895153_A"
    args$ref_db <- glue("{here::here()}/001_data_local/interactions_db_v2/ref_db.rds")
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
pacman::p_load(reshape2)

log_info("Load reference database...")
ref_db <- readRDS(args$ref_db) %>%
    select(interaction, simple_interaction, complex_interaction)

log_info("Load Cell2Cell output...")
interactions <- data.frame(fread(args$input_interactions, header = TRUE, sep = "\t"),
    check.names = FALSE
)
log_info("Formatting results...")
interactions <- interactions %>%
    rowwise() %>%
    mutate(interaction = V1 %>% str_replace_all("\\(|,|\\'|\\)", "") %>%
        str_replace_all(" ", "_")) %>%
    ungroup() %>%
    select(-V1) %>%
    reshape2::melt(id.vars = c("interaction"), variable.name = "source_target", value.name = "pval") %>%
    mutate(source_target = str_replace_all(source_target, ";", "__"), method = "Cell2Cell", Sample = args$sample_id) %>%
    left_join(ref_db, by = "interaction") %>%
    select(-interaction)


log_info("Save output...")
saveRDS(interactions,
    file = glue("{args$output_dir}/cell2cell__{args$sample}__postproc.rds")
)

log_info("COMPLETED!")
