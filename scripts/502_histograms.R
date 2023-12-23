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
        description = "Ligand/receptor expression vs background",
    )
    parser$add_argument("--gene_expression", type = "character", help = "Path to gene expression matrix")
    parser$add_argument("--is_stringent", type = "numeric", help = "Stringency of interactions", default = 1)
    parser$add_argument("--interactions", type = "character", help = "Path to interactions (post-filtering)")
    parser$add_argument("--binwidth", type = "numeric", help = "Bin width", default = 0.3)
    parser$add_argument("--max_q", type = "numeric", help = "Max. quantile", default = 0.85)
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$is_stringent <- 1
    args$gene_expression <- glue("{here::here()}/final_output/neuron_tumor/002_prep_sample//6419_enhancing_border__prepped.rds")
    args$interactions <- glue("{here::here()}/final_output/macrophage_tumor/402_post_filtering/ccis_post_filtering__stringency_{args$is_stringent}_groupby_1.rds")
    args$output_dir <- glue("{here::here()}/final_output/macrophage_tumor")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- paste0(args$output_dir, "/502_histograms")
create_dir(output_dir)

# Load additional libraries
pacman::p_load(ggplot2, Seurat, ggrepel, gridExtra, grid, ggforce)

log_info("Load gene expression...")
obj <- readRDS(args$gene_expression)

rename_dict <- c("AC-like" = "Malignant_AC", "MES-like" = "Malignant_MES", "NPC & OPC-like" = "Malignant_NPC_&_OPC", "TAM-BDM" = "Macrophage")


# TODO Should be replaced with just one metadata object
metadata_parse <- readRDS(glue("{here::here()}/output_tmp/data/metadata__parsebio.rds")) %>% mutate(CellClass_L4 = NA, Dataset = "ParseBio")
metadata_10x <- readRDS(glue("{here::here()}/output_tmp/data/metadata__10x.rds")) %>%
    select(-c(SampleID)) %>%
    mutate(Sample = str_replace_all(Sample, "__", "_"), Dataset = "Multiome 10x")

all_meta <- rbind(metadata_parse, metadata_10x) %>%
    select(Sample, Batch, Patient, Region, Dataset) %>%
    distinct()


log_info("Extract sample information...")
current_sample_id <- str_remove(get_name(args$gene_expression), "__prepped")
current_region <- all_meta %>%
    filter(Sample == current_sample_id) %>%
    pull(Region) %>%
    unique()

log_info(glue("\n\nSample ID: {current_sample_id}\nRegion: {current_region}\n\n"))


log_info("Load reference database of interactions...")
db_ref <- readRDS(glue("{here::here()}/data/interactions_db/interactions_ref.rds"))

log_info("Load interactions...")
all_interactions <- readRDS(args$interactions) %>% separate(interaction, into = c("ligand", "receptor"), sep = " - ", remove = FALSE)

log_info("Extract source-targets of interest...")
source_targets_oi <- all_interactions %>%
    ungroup() %>%
    select(source_target) %>%
    separate(source_target, into = c("source", "target"), sep = "__", remove = FALSE) %>%
    unique()

avg_expr_by_group <- as.data.frame(AverageExpression(obj, assay = "RNA", group.by = "custom_annot")$RNA)

colnames(avg_expr_by_group) <- str_replace_all(colnames(avg_expr_by_group), rename_dict)

for (i in seq_len(nrow(source_targets_oi))) {
    current_source <- source_targets_oi[i, ]$source
    current_target <- source_targets_oi[i, ]$target
    log_info(glue("Processing {i} of {nrow(source_targets_oi)}: {current_source}-{current_target}..."))

    log_info("Filtering interactions...")
    interactions <- all_interactions %>% filter(Region == current_region, source_target == glue("{current_source}__{current_target}"))

    log_info("Extract ligands + receptors from interactions...")
    ligands_oi <- interactions %>%
        filter(source == current_source) %>%
        pull(ligand) %>%
        unique()
    receptors_oi <- interactions %>%
        filter(target == current_target) %>%
        pull(receptor) %>%
        unique()

    log_info("Average expression for source (ligands) and targets (receptors) separately...")
    obj_source_avg_expr <- data.frame(expr = avg_expr_by_group[, current_source]) %>%
        mutate(gene = rownames(avg_expr_by_group)) %>%
        filter(expr > 0)

    obj_target_avg_expr <- data.frame(expr = avg_expr_by_group[, current_target]) %>%
        mutate(gene = rownames(avg_expr_by_group)) %>%
        filter(expr > 0)

    log_info("Select ligands + receptors of interest from interactions for labelling...")
    source_annot <- obj_source_avg_expr %>%
        filter(gene %in% ligands_oi) %>%
        mutate(y = 0)
    target_annot <- obj_target_avg_expr %>%
        filter(gene %in% receptors_oi) %>%
        mutate(y = 0)

    log_info("Visualize...")
    source_plots <- plot_expr_hist_ccis(obj_source_avg_expr, source_annot, args)
    target_plots <- plot_expr_hist_ccis(obj_target_avg_expr, target_annot, args)
    p_combi_log10 <- grid.arrange(
        source_plots[[2]] + labs(title = glue("Source: {current_source}")), target_plots[[2]] + labs(title = glue("Target: {current_target}")),
        ncol = 2,
        bottom = textGrob("Average expression", gp = gpar(fontsize = 14, fontface = "bold")),
        left = textGrob("log10(Number of cells)", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
    )

    p_combi <- grid.arrange(
        source_plots[[1]] + labs(title = glue("Source: {current_source}")), target_plots[[1]] + labs(title = glue("Target: {current_target}")),
        ncol = 2,
        bottom = textGrob("Average expression", gp = gpar(fontsize = 14, fontface = "bold")),
        left = textGrob("Number of cells", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
    )
    log_info("Save plot...")
    ggsave(plot = p_combi_log10, file = glue("{output_dir}/{current_sample_id}__{current_source}_{current_target}_hist_log10__stringency_{args$is_stringent}.pdf"), width = 15, height = 8, dpi = 300)

    ggsave(plot = p_combi, file = glue("{output_dir}/{current_sample_id}__{current_source}_{current_target}_hist__stringency_{args$is_stringent}.pdf"), width = 15, height = 8, dpi = 300)
}

log_info("COMPLETED!")
