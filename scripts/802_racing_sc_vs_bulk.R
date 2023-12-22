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
    args$min_samples <- 3
    args$output_dir <- glue("{here::here()}/output/RaCInG/figures")
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
pacman::p_load_gh("jokergoo/ComplexHeatmap", "jmw86069/colorjam")
pacman::p_load(VennDiagram)


racing_tam_tumor <- readRDS(glue("{here::here()}/output/RaCInG/top_interactions__tam_tumor.rds"))
racing_tumor_tam <- readRDS(glue("{here::here()}/output/RaCInG/top_interactions__tumor_tam.rds"))

racing_tumor_tam <- racing_tumor_tam %>% mutate(source_target = paste0(top_interactions, "__TAM-BDM"))
racing_tam_tumor <- racing_tam_tumor %>% mutate(source_target = paste0("TAM-BDM__", top_interactions))

racing <- rbind(racing_tumor_tam, racing_tam_tumor) %>%
    mutate(method = "bulk") %>%
    select(interaction, source_target, method)
# ---- PREPARE SC-RNASEQ RESULTS ---- #
# Contains both 10x and parsebio data
db <- readRDS(glue("{here::here()}/output/RaCInG/all_weights__tam_tumor.rds")) %>% unite("interaction", c("ligand", "receptor"), sep = "__")

interactions_ref <- readRDS(glue("{here::here()}/data/interactions_db/interactions_ref.rds"))

pairs_oi <- read.table(glue("{here::here()}/data/source_target_oi.txt")) %>% pull(V1)

lit_obj <- readRDS(glue("{here::here()}/data/literature_findings_v3.rds")) %>% rename(setname = source_target)

common_db <- str_replace_all(intersect(db$interaction, interactions_ref$interaction), "__", " - ")

ccis_combined <- readRDS(glue("{here::here()}/output/tam-bdm_tumor/min_200_cells/all_ccis_combined_filtered.rds")) %>%
    filter(interaction %in% common_db, source_target %in% pairs_oi)

ccis_oi <- ccis_combined %>%
    # Binarize
    filter(n >= args$min_samples) %>%
    mutate(method = "sc") %>%
    ungroup() %>%
    select(interaction, source_target, method)

racing <- racing %>% filter(interaction %in% common_db)

combined <- rbind(ccis_oi, racing) %>% mutate(setname = case_when(source_target %in% pairs_oi[1:3] ~ "TAM-BDM - Malignant", source_target %in% pairs_oi[4:6] ~ "Malignant - TAM-BDM"))

# ---- PLOT ---- #
for (pair in c("TAM-BDM - Malignant", "Malignant - TAM-BDM")) {
    # Per pair set
    sc <- combined %>%
        filter(setname == pair, method == "sc") %>%
        pull(interaction) %>%
        unique()
    bulk <- combined %>%
        filter(setname == pair, method == "bulk") %>%
        pull(interaction) %>%
        unique()
    data <- list(sc = sc, bulk = bulk)
    create_venndiagram(x = data, filename = glue("{args$output_dir}/sc_bulk_venn__{str_replace_all(pair, c(' - ' = '_'))}.png"), category.names = names(data), main = pair)

    # Intersection
    intersection <- intersect(sc, bulk)
    log_info(glue("Number of interactions in intersection: {length(intersection)}"))
    log_info(glue("Common interactions: {paste(intersection, collapse = ', ')}"))
}

for (pair in pairs_oi) {
    # Per pair set
    sc <- combined %>%
        filter(source_target == pair, method == "sc") %>%
        pull(interaction) %>%
        unique()
    bulk <- combined %>%
        filter(source_target == pair, method == "bulk") %>%
        pull(interaction) %>%
        unique()
    data <- list(sc = sc, bulk = bulk)
    create_venndiagram(x = data, filename = glue("{args$output_dir}/sc_bulk_venn__{str_replace_all(pair, c(' - ' = '_'))}.png"), category.names = names(data), main = pair)

    # Intersection
    log_info(glue("Pair: {pair}"))
    intersection <- intersect(sc, bulk)
    log_info(glue("Number of interactions in intersection: {length(intersection)}"))
    log_info(glue("Common interactions: {paste(intersection, collapse = ', ')}"))
    print(" ")
}

combined <- combined %>%
    left_join(lit_obj %>% distinct(setname, interaction) %>% mutate(found_previous = 1), by = c("setname", "interaction")) %>%
    mutate(is_identified = 1) %>%
    pivot_wider(names_from = method, values_from = is_identified)

write.table(combined, file = glue("{here::here()}/output/RaCInG/comparison_top25bulk_vs_sc.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
