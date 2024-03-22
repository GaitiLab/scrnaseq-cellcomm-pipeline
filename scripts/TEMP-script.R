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
    args$output_dir <- glue("{here::here()}/output/")
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


obj <- readRDS(glue("{here::here()}/data/interactions_db/ref_db.rds")) %>%
    rowwise() %>%
    mutate(
        source_genesymbol_ordered = paste0(sort(str_split(source_genesymbol, "_", simplify = TRUE)), collapse = "_"),
        target_genesymbol_ordered = paste0(sort(str_split(target_genesymbol, "_", simplify = TRUE)), collapse = "_")
    )

obj %>% filter(str_detect(source_genesymbol, "NLGN"))

obj %>% nrow()

# obj %>%
#     distinct(complex_interaction) %>%
#     nrow()

# pathways_cellchat_og <- read.csv(glue("{here::here()}/000_misc_local/cellchat_pathway.csv"))


# pathways_cellchat <- pathways_cellchat_og %>%
#     mutate(
#         source_genesymbol = str_replace_all(ligand.symbol, ", ", "_"),
#         target_genesymbol = str_replace_all(receptor.symbol, ", ", "_")
#     ) %>%
#     mutate(
#         source_genesymbol = ifelse(source_genesymbol == "", toupper(ligand), source_genesymbol),
#         target_genesymbol = ifelse(target_genesymbol == "", toupper(receptor), target_genesymbol)
#     ) %>%
#     rowwise() %>%
#     mutate(
#         source_genesymbol_ordered = paste0(sort(str_split(source_genesymbol, "_", simplify = TRUE)), collapse = "_"),
#         target_genesymbol_ordered = paste0(sort(str_split(target_genesymbol, "_", simplify = TRUE)), collapse = "_")
#     ) %>%
#     select(-X, -interaction_name, -ligand, -receptor, -interaction_name_2, -source_genesymbol, -target_genesymbol)



# obj_merged <- obj %>% left_join(pathways_cellchat, by = c("source_genesymbol_ordered", "target_genesymbol_ordered"))

# obj_merged %>% nrow()
# obj_merged %>%
#     distinct(complex_interaction) %>%
#     nrow()


# saveRDS(obj_merged, glue("{here::here()}/data/interactions_db/ref_db_with_pathway.rds"))
