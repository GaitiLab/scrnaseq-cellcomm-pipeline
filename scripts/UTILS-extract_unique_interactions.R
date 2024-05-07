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
        description = "Extract unique interactions",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("/Users/joankant/Library/CloudStorage/OneDrive-UHN/Spatial_GBM/Analysis/WIP/scRNAseq/CCI")
    args$input_file <- "output/CCI_CellClass_L2_2_reassigned_samples_confident_only/402_aggregation/402c_aggregation_integration.rds"
    args$condition_varname <- "Region"
    args$output_name <- "CCI_CellClass_L2_2_reassigned_samples_confident_only_unique_interactions_neuron_invasive_high"
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
pacman::p_load(xlsx)

output_filename <- glue("{args$output_dir}/{args$output_name}.xlsx")
if (file.exists(output_filename)) {
    file.remove(output_filename)
}

log_info("Load data + formatting...")
obj <- readRDS(args$input_file) %>%
    # Additional filtering on p-value (aggregate rank)
    filter(pval < 0.05) %>%
    separate(source_target, c("source", "target"), sep = "__", remove = FALSE)


obj_undirected <- obj %>%
    # Remove direction by sorting source-target alphabetically
    rowwise() %>%
    mutate(source_target_undirected = paste0(sort(c(source, target)), collapse = "__")) %>%
    # Remove duplicate interactions, when they are found in both directions only keep 1
    distinct(!!sym(args$condition_varname), source_target_undirected, complex_interaction, .keep_all = FALSE)

celltype_pairs_oi <- c(
    "Neuron__Invasive-high OPC/NPC1",
    "Invasive-high OPC/NPC1__Neuron",
    "Neuron__Progenitor_like",
    "Progenitor_like__Neuron"
)

obj_filtered <- obj_undirected %>%
    filter(!!sym(args$condition_varname) == "PT", source_target_undirected %in% celltype_pairs_oi) %>%
    mutate(dummy = 1) %>%
    pivot_wider(names_from = source_target_undirected, values_from = dummy, values_fill = 0) %>%
    # Mark unique interactions
    mutate(is_unique = ifelse((`Invasive-high OPC/NPC1__Neuron` == 1) & (Neuron__Progenitor_like == 0), 1, 0))


unique_interactions <- obj_filtered %>%
    filter(is_unique == 1) %>%
    pull(complex_interaction)


unique_interactions_df <- obj %>% filter(complex_interaction %in% unique_interactions, source_target %in% celltype_pairs_oi[1:2], !!sym(args$condition_varname) == "PT")


log_info("Write data to Excel...")
write.xlsx(unique_interactions_df, output_filename,
    col.names = TRUE, row.names = FALSE, append = TRUE
)
log_info("COMPLETED!")
