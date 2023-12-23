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
    args$meta_parse <- glue("{here::here()}/output/parsebio_new/000_misc/metadata.rds")
    args$meta_10x <- glue("
    {here::here()}/output/gbm10x_new/000_misc/metadata.rds")
    args$output_dir <- glue("{here::here()}/output/HE_images")
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
pacman::p_load(readxl, ggplot2, ggpubr)

gbm360 <- read_excel("/Users/joankant/Library/CloudStorage/OneDrive-UHN/GBM/pathology_images_checklist.xlsx", sheet = "cell fraction")
head(gbm360)

gbm360 %>% pull(cell_type) %>% unique()

meta_10x <- readRDS(args$meta_10x)
gbm10x_cell_types <- paste0(meta_10x %>% pull(Malignant.subtyping.L2) %>% unique()
, collapse = ', ')
log_info(glue("Cell types: {gbm10x_cell_types}"))

# ---- Deal with 10x data ---- #
gbm10x_annot <- meta_10x %>% mutate(cell_type = case_when(
    Malignant.subtyping.L2 == "OPC" ~ "Normal",
    Malignant.subtyping.L2 == "Malignant.AC" ~ "AC",
    Malignant.subtyping.L2 == "Malignant.NPC1" | Malignant.subtyping.L2 == "Malignant.NPC2" ~ "NPC",
    Malignant.subtyping.L2 == "Malignant.MES_ASTROCYTE"  | Malignant.subtyping.L2 == "Malignant.MES_INT" ~ "MESlike",
    Malignant.subtyping.L2 == "Malignant.MES_HYPOXIA" ~ "MEShypoxia",
    Malignant.subtyping.L2 == "Malignant.OPC" ~ "OPC",
    Malignant.subtyping.L2 == "NoSignals" ~ NA,
    TRUE ~ "Normal")
) %>% filter(!is.na(cell_type)) %>% select(cell_type, sample)
n_by_cell_type <- gbm10x_annot %>% group_by(sample, cell_type) %>% reframe(n_cells= n())

gbm10x_fracs <- n_by_cell_type %>% group_by(sample) %>% mutate(prop = n_cells/ sum(n_cells))

# ---- Merge data ---- #
sc_fracs <- gbm10x_fracs %>% mutate(Percentage = prop * 100, method = "scRNA-seq") %>% select(sample, cell_type, Percentage, method)

gbm360 <- gbm360 %>% mutate(method = "gbm360") %>% select(sample, cell_type, Percentage, method)

common_samples <- gbm360 %>% pull(sample) %>% intersect(sc_fracs %>% pull(sample))

merged <- rbind(gbm360, sc_fracs) %>% filter(sample %in% common_samples)

head(merged)

p <- ggplot(merged, aes(fill=cell_type, y=Percentage, x=method)) +
    geom_bar(position="stack", stat="identity") +
    facet_wrap(~sample)
ggsave(p, file = glue("{args$output_dir}/gbm360_vs_scrnaseq.pdf"), width = 10, height = 10, dpi = 300)
