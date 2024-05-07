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
        description = "Correlation analysis bulkRNAseq Visualization",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/corr_analysis_bulk_rnaseq_unique_interactions")
    args$prefix <- "TCGA"
    args$suffix <- "CIBERSORT"
    # args$input_file <- "output/corr_analysis_bulk_rnaseq/GLASS_paired_corr_boot.rds"
    args$input_file <- glue("output/corr_analysis_bulk_rnaseq_unique_interactions/{args$prefix}_corr_boot_{args$suffix}.rds")
    args$interactions <- "000_misc_local/unique_interactions.xlsx"
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
pacman::p_load(ggplot2, readxl)
obj <- readRDS(args$input_file)

# Interactions
interactions_subset <- read_excel(args$interactions) %>% arrange(pval) %>% 
    separate(complex_interaction, into = c("ligand_complex", "receptor_complex"), sep = "__")

neuron_ligands <- interactions_subset %>% filter(startsWith(source_target, "Neuron__")) %>% pull(ligand_complex) %>% str_split(., ":") %>% unlist() %>% unique()
neuron_receptors <-  interactions_subset %>% filter(startsWith(source_target, "__Neuron")) %>% pull(ligand_complex) %>% str_split(., ":") %>% unlist() %>% unique()

invasive_ligands <- interactions_subset %>% filter(startsWith(source_target, "Neuronal OPC-like")) %>% pull(ligand_complex) %>% str_split(., ":") %>% unlist() %>% unique()
invasive_receptors <-  interactions_subset %>% filter(startsWith(source_target, "Neuronal OPC-like")) %>% pull(ligand_complex) %>% str_split(., ":") %>% unlist() %>% unique()

invasive_genes <- unique(c(invasive_ligands, invasive_receptors))

color_palette <- load_color_palette("ATAC")
names(color_palette) <- str_replace_all(names(color_palette), "Neuronal.OPC.like", "Invasive-high OPC/NPC1")


if (str_detect(args$suffix, "CIBERSORT")) {
    scalop_res <- obj$scalop %>%
        mutate(pair = str_replace_all(pair, "Neuronal.OPC.like", "Invasive-high OPC/NPC1"), gene = factor(gene)) %>%
        separate(pair, into = c("state", "dummy"), sep = "__", remove = FALSE) %>%
        mutate(state = factor(state, names(color_palette))) %>%
        rowwise() %>%
        filter(str_detect(state, "Malignant") | state == "Invasive-high OPC/NPC1")
    ssgsea_res <- obj$ssgsea %>%
        mutate(pair = str_replace_all(pair, "Neuronal.OPC.like", "Invasive-high OPC/NPC1"), gene = factor(gene)) %>%
        separate(pair, into = c("state", "dummy"), sep = "__", remove = FALSE) %>%
        rowwise() %>%
        filter(str_detect(state, "Malignant") | state == "Invasive-high OPC/NPC1") %>%
        mutate(state = factor(state, names(color_palette)))
} else {
    scalop_res <- obj$scalop %>%
        mutate(gene = factor(gene), state = factor(state))
    ssgsea_res <- obj$ssgsea %>%
        mutate(gene = factor(gene), state = factor(state))
}

log_info("Visualize...")
p_scalop <- ggplot(data = scalop_res, aes(x = gene, y = corr_mean, fill = state)) +
    geom_bar(
        stat = "identity", color = "black",
        position = "dodge"
    ) +
    guides(fill = guide_legend(override.aes = list(size = 7), title = "Label")) +
    geom_errorbar(aes(ymin = corr_mean - corr_se, ymax = corr_mean + corr_se),
        width = .2,
        position = position_dodge(.9)
    ) +
    scale_fill_manual(values = color_palette) +
    GBM_theme() +
    labs(x = "Receptor", y = "Correlation") +
    theme(axis.text.x = element_blank())
p_scalop

p_ssgsea <- ggplot(data = ssgsea_res, aes(x = gene, y = corr_mean, fill = state)) +
    geom_bar(
        stat = "identity", color = "black",
        position = position_dodge()
    ) +
    guides(fill = guide_legend(override.aes = list(size = 7), title = "Label")) +
    geom_errorbar(aes(ymin = corr_mean - corr_se, ymax = corr_mean + corr_se),
        width = .2,
        position = position_dodge(.9)
    ) +
    GBM_theme() +
    labs(x = "Receptor", y = "Correlation")
p_ssgsea


if (str_detect(args$prefix, "GLASS")) {
    p_ssgsea <- p_ssgsea + ggh4x::facet_wrap2(gene ~ Group, scales = "free_x")
    p_scalop <- p_scalop + ggh4x::facet_wrap2(gene ~ Group, scales = "free_x")
} else {
    p_ssgsea <- p_ssgsea + ggh4x::facet_wrap2(. ~ gene, ncol = 8, scales = "free_x")
    p_scalop <- p_scalop + ggh4x::facet_wrap2(. ~ gene, ncol = 8, scales = "free_x")
}

log_info("Save plot...")

ggsave(p_scalop, filename = glue("{args$output_dir}/{get_name(args$input_file)}_scalop.pdf"), width = 45, height = 30)

ggsave(p_ssgsea, filename = glue("{args$output_dir}/{get_name(args$input_file)}_ssgsea.pdf"), width = 45, height = 30)
