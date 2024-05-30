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
    args$output_dir <- glue("{here::here()}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/302_postproc_cell2cell")
    args$input_interactions <- glue("{here::here()}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/202_cci_cell2cell/cell2cell__6509_cortex.csv")
    args$ref_db <- glue("{here::here()}/data/interactions_db/ref_db.rds")
    args$sample_id <- "6509_cortex"
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

# Setting sample + run id depending on downsampling or not
split_sample_id <- str_split(args$sample_id, "__", simplify = TRUE)
sample_id <- split_sample_id[1]
run_id <- NA
if (length(split_sample_id) > 2) {
    run_id <- split_sample_id[3]
}

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
    mutate(source_target = str_replace_all(source_target, ";", "__"), method = "Cell2Cell", Sample = sample_id, run_id = run_id) %>%
    left_join(ref_db, by = "interaction") %>%
    select(-interaction)

#       source_target        pval    method         Sample run_id simple_interaction  complex_interaction
# 1 Neuron__Microglia 0.461141321 Cell2Cell 6234_2895153_A      1        NRG2__ERBB2    NRG2__ERBB2:ERBB3
# 2 Neuron__Microglia 0.461141321 Cell2Cell 6234_2895153_A      1        NRG2__ERBB2    NRG2__ERBB2:ERBB4
# 3 Neuron__Microglia 0.999083070 Cell2Cell 6234_2895153_A      1         EREG__EGFR     EREG__EGFR:ERBB2
# 4 Neuron__Microglia 0.999083070 Cell2Cell 6234_2895153_A      1        EREG__ERBB2    EREG__ERBB2:ERBB4
# 5 Neuron__Microglia 0.177891249 Cell2Cell 6234_2895153_A      1      SEMA7A__ITGA1  SEMA7A__ITGA1:ITGB1
# 6 Neuron__Microglia 0.008461169 Cell2Cell 6234_2895153_A      1       NFASC__CNTN1 NFASC__CNTN1:CNTNAP1

log_info("Save output...")
saveRDS(interactions,
    file = glue("{args$output_dir}/cell2cell__{args$sample}__postproc.rds")
)

log_info("COMPLETED!")
