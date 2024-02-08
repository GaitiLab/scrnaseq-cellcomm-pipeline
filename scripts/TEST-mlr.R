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
        description = "Test relationship with MLR",
    )
    parser$add_argument("--annot", type = "character", help = "Annotation to use")
    parser$add_argument("--type_of_voting", type = "character", help = "'lenient_voting' or 'stringent_voting'")
    parser$add_argument("--meta", type = "character", help = "metadata object RDS")
    parser$add_argument("--interactions_n_cells", type = "character", help = "Number of interactions per cells (output from TEST-ncells_impact.R)")
    parser$add_argument("--potential_interactions", type = "character", help = "Number of (detected/potential) interactions per sample (output from TEST-potential_vs_detected)")
    parser$add_argument("--genes_per_class", type = "character", help = "Number of genes expressed per class (output from TEST-potential_vs_detected)")
    parser$add_argument("--interactions_db", type = "character", help = "Database of interactions")


    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L1"
    args$type_of_voting <- "lenient_voting"
    args$output_dir <- glue("{here::here()}/output/{args$annot}/TEST-mlr")
    args$meta <- glue("{here::here()}/output/CCI_CellClass_L2/000_data/bw_gbm_regional_study__metadata.rds")
    args$interactions_n_cells <- glue("{here::here()}/output/{args$annot}/TEST-ncells_impact/all_n_interactions_by_sample_{args$type_of_voting}.rds")
    args$potential_interactions <- glue("{here::here()}/output/{args$annot}/TEST-potential_vs_detected/all_n_interactions_by_sample_{args$type_of_voting}.rds")
    args$genes_per_class <- glue("{here::here()}/output/{args$annot}/TEST-potential_vs_detected/genes_per_class_{args$type_of_voting}.rds")
    args$interactions_db <- glue("{here::here()}/data/interactions_db/ref_db.rds")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- glue("{args$output_dir}/{args$type_of_voting}")
create_dir(output_dir)
set.seed(123)
# Load additional libraries
pacman::p_load(ggplot2, ggcorrplot, ggpubr, ggtext)

meta <- readRDS(args$meta)
db <- readRDS(args$interactions_db)

db_ligands <- str_split(db %>%
    pull(source_genesymbol) %>%
    unique(), "_", simplify = TRUE)[, 1] %>% unique()
db_receptors <- str_split(db %>%
    pull(target_genesymbol) %>%
    unique(), "_", simplify = TRUE)[, 1] %>% unique()

# Related to expression / number of cells
interactions_n_cells <- readRDS(args$interactions_n_cells)
potential_interactions <- readRDS(args$potential_interactions)

# Expressed ligands/receptors per class
genes_per_class <- readRDS(args$genes_per_class)
genes_per_row <- genes_per_class %>%
    pull(expressed_genes) %>%
    unlist() %>%
    str_split(., ", ")

# Determine number of ligands/receptors
n_lr <- do.call(rbind, lapply(genes_per_row, function(x) {
    n_ligands <- length(intersect(x, db_ligands))
    n_receptors <- length(intersect(x, db_receptors))
    return(c(n_ligands, n_receptors))
}))
n_lr <- data.frame(n_lr)
colnames(n_lr) <- c("n_ligands", "n_receptors")
genes_per_class <- cbind(genes_per_class, n_lr) %>% rename(cell_type = id)
genes_per_class_subset <- genes_per_class %>%
    select(cell_type, Sample, Region_Grouped, n_ligands, n_receptors)

combi <- merge(interactions_n_cells, potential_interactions)
out <- combi %>%
    left_join(genes_per_class_subset, by = c("source" = "cell_type", "Sample", "Region_Grouped")) %>%
    rename(n_ligands_source = n_ligands, n_receptors_source = n_receptors) %>%
    left_join(genes_per_class_subset, by = c("target" = "cell_type", "Sample", "Region_Grouped")) %>%
    rename(n_ligands_target = n_ligands, n_receptors_target = n_receptors) %>%
    mutate(n_lr_source = n_ligands_source + n_receptors_source, n_lr_target = n_ligands_target + n_receptors_target) %>%
    mutate(Region_Grouped = as.numeric(as.factor(Region_Grouped)))

predictors <- out %>% select(Region_Grouped, n_interactions, n_lr_source, n_lr_target, n_cells_source, n_cells_target, n_possible_interactions, n_ligands_source, n_receptors_source, n_ligands_target, n_receptors_target)
predictors_long <- predictors %>% reshape2::melt()

# Attempt 1
model_1 <- lm(n_interactions ~ Region_Grouped + n_lr_source + n_lr_target + n_cells_source + n_cells_target + n_ligands_source + n_receptors_source + n_ligands_target + n_receptors_target + n_possible_interactions, data = out)
model_1_residuals <- data.frame(model_1$residuals)

plt_residuals <- ggplot(model_1_residuals) +
    geom_histogram(aes(model_1.residuals)) +
    custom_theme() +
    labs(x = "Residuals", y = "Count")

plt_normality <- ggplot(data = model_1_residuals, aes(sample = model_1.residuals)) +
    stat_qq() +
    stat_qq_line() +
    custom_theme() +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

plt_quality <- ggarrange(plt_residuals, plt_normality)
ggsave(plot = plt_quality, filename = glue("{output_dir}/plt_model_1.pdf"), width = 15, height = 5)

sink(glue("{output_dir}/model1_all_predictors.txt"))
print(summary(model_1))

# Distribution
ggplot(data = predictors_long, aes(value)) +
    geom_histogram() +
    facet_wrap(~variable, scales = "free") +
    custom_theme()

# Multi-collinearity
corr_mat <- round(predictors %>%
    cor(use = "pairwise.complete.obs", method = "spearman"), 2)

plt_corr <- ggcorrplot(corr_mat,
    hc.order = TRUE, type = "lower",
    lab = TRUE
) + labs(x = "", y = "") + guides(fill = guide_colourbar(title = "Spearman Coefficient", nbin = 3, barwidth = 10, barheight = 1)) + custom_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2))
plt_corr
ggsave(plot = plt_corr, filename = glue("{output_dir}/corr_mat.pdf"), width = 8, height = 8)

# Removing highly correlated predictors (features)
predictors_subset <- predictors %>% select(Region_Grouped, n_cells_source, n_cells_target, n_ligands_source, n_ligands_target, n_possible_interactions)

# Attempt 2
model_2 <- lm(n_interactions ~ Region_Grouped + n_cells_source + n_cells_target + n_ligands_source + n_ligands_target + n_possible_interactions, data = out)
model_2_residuals <- data.frame(model_2$residuals)

plt_residuals <- ggplot(model_2_residuals) +
    geom_histogram(aes(model_2.residuals)) +
    custom_theme() +
    labs(x = "Residuals", y = "Count")

plt_normality <- ggplot(data = model_2_residuals, aes(sample = model_2.residuals)) +
    stat_qq() +
    stat_qq_line() +
    custom_theme() +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

plt_quality <- ggarrange(plt_residuals, plt_normality)
ggsave(plot = plt_quality, filename = glue("{output_dir}/plt_model_2.pdf"), width = 15, height = 5)

sink(glue("{output_dir}/model2_remove_highly_correlated.txt"))
print(summary(model_2))


# Anova test
comp <- anova(model_1, model_2)
print(comp$`Pr(>F)`)[2]
# p > 0.05 -> model is an improvement

# Attempt 3
model_3 <- lm(n_interactions ~ Region_Grouped + n_cells_source + n_cells_target + n_possible_interactions, data = out)
model_3_residuals <- data.frame(model_3$residuals)

plt_residuals <- ggplot(model_3_residuals) +
    geom_histogram(aes(model_3.residuals)) +
    custom_theme() +
    labs(x = "Residuals", y = "Count")

plt_normality <- ggplot(data = model_3_residuals, aes(sample = model_3.residuals)) +
    stat_qq() +
    stat_qq_line() +
    custom_theme() +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

plt_quality <- ggarrange(plt_residuals, plt_normality)
ggsave(plot = plt_quality, filename = glue("{output_dir}/plt_model_3.pdf"), width = 15, height = 5)

comp <- anova(model_2, model_3)
print(comp$`Pr(>F)`)[2]

sink(glue("{output_dir}/model3_remove_insignif_pred.txt"))
print(summary(model_3))
