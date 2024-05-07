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
        description = "Score pseudobulking",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/002_WIP_FIGURES_local/pseudobulk_scRNAseq_subtyping_Verhaak")
    # args$input_file <- "output/pseudobulk_scoring/pseudobulking_scored_ssgsea_w_filtered_signatures.rds"

    # args$input_file <- "output/pseudobulk_scoring/pseudobulk_scored_scalop_w_filtered_signatures.rds"
    # args$input_file <- "output/pseudobulk_scoring/pseudobulk_scored_scalop_w_full_signatures.rds"

    args$method <- "scalop"
    args$suffix <- "all_gmt"

    args$input_file <- glue("output/pseudobulk_scoring/pseudobulking_scored_{args$method}_{args$suffix}_signatures.rds")
    args$label <- "VERHAAK"

    # args$input_file <- glue("output/pseudobulk_scoring/pseudobulking_scored_{args$method}_public_signatures.rds")
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
pacman::p_load(GSVA, readxl, boot, scalop, Hmisc, corrplot)
pacman::p_load(ggpubr)

obj <- readRDS(args$input_file)
head(obj)
if (args$method == "scalop") {
    obj <- t(obj)
}


sample_meta <- readRDS("000_misc_local/gbm_regional_study__metadata.rds") %>%
    remove_rownames() %>%
    select(Sample, Platform, Region) %>%
    distinct()

obj_long <- obj %>%
    as.data.frame() %>%
    rownames_to_column("Label") %>%
    # Rename our signature to updated label
    mutate(Label = ifelse(Label == "Neuronal.OPC.like", "Invasive-high OPC/NPC1", Label)) %>%
    reshape2::melt(id.vars = "Label", value.name = "score", variable.name = "Sample")

sample_subtype <- obj_long %>%
    group_by(Sample) %>%
    top_n(1, score) %>%
    rename(subtype = Label, subtype_score = score) %>%
    left_join(sample_meta)

merged <- obj_long %>%
    left_join(sample_meta) %>%
    filter(str_detect(Label, args$label)) %>%
    mutate(
        Region = factor(Region, levels = names(GBMutils::load_color_palette("Region"))),
        Label = factor(str_to_title(str_remove(Label, "VERHAAK_GLIOBLASTOMA_")), levels = names(GBMutils::load_color_palette("Verhaak")))
    )
comparisons <- generate_pairs(names(GBMutils::load_color_palette("Region"))) %>% str_split(., "__")

unique_labels <- merged %>%
    pull(Label) %>%
    unique()
for (label in unique_labels) {
    p <- ggboxplot(merged %>% filter(Label == label), x = "Region", y = "score", color = "Region") + stat_compare_means(comparisons = comparisons) +
        geom_point(data = merged %>% filter(Label == label), aes(x = Region, y = score, color = Region), show.legend = FALSE, position = position_jitter(width = 0.1, height = 0.1)) +
        GBM_theme() +
        labs(x = "Region", y = "Gene signature score", title = label, subtitle = glue("n={merged %>% pull(Sample) %>% unique() %>% length()}")) +
        scale_color_manual(
            values = GBMutils::load_color_palette("Region"),
            guide = guide_legend(override.aes = list(size = 12))
        ) +
        scale_y_continuous(breaks = scales::pretty_breaks())
    p
    label <- str_replace(label, "/", "_")
    ggsave(p, filename = glue("{args$output_dir}/boxplots_{label}__{args$method}__{args$label}_signatures.pdf"), width = 10, height = 7, dpi = 300)
}

# Distribution of assigned subtypes across regions
subtypes <- sample_subtype %>%
    mutate(Region = factor(Region, levels = names(GBMutils::load_color_palette("Region"))), subtype = factor(str_to_title(str_remove(subtype, "VERHAAK_GLIOBLASTOMA_")), levels = names(GBMutils::load_color_palette("Verhaak"))))
saveRDS(subtypes, file = glue("{args$output_dir}/distribution__{args$label}__{args$method}__{args$suffix}_signatures.rds"))

subtypes_by_region <- sample_subtype %>%
    group_by(Region) %>%
    count(subtype)

# Statistical test
subtypes_by_region_wide <- subtypes_by_region %>%
    pivot_wider(names_from = "Region", values_from = "n", values_fill = 0) %>%
    column_to_rownames("subtype")

# Chisq test: are row and column vars independent.
prop.table(table(subtypes$subtype, subtypes$Region), margin = 2)
chisq.test(subtypes$subtype, subtypes$Region)
# Result:
#         Pearson's Chi-squared test

# data:  subtypes$subtype and subtypes$Region
# X-squared = 25.329, df = 6, p-value = 0.0002968

# Warning message:
# In chisq.test(subtypes$subtype, subtypes$Region) :
#   Chi-squared approximation may be incorrect
# insufficient number of counts, therefore results not reliable

# Fisher test
stat_test <- stats::fisher.test(subtypes_by_region_wide)
print(stat_test)

#         Fisher's Exact Test for Count Data

# data:  subtypes_by_region_wide
# p-value = 8.386e-05
# alternative hypothesis: two.sided

pval <- formatC(stat_test$p.value, format = "e", digits = 2)

p <- ggplot(data = subtypes_by_region, aes(x = Region, y = n, fill = subtype)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(x = "Region", y = "Number of samples", subtitle = glue("n={sample_subtype %>% pull(Sample) %>% unique() %>% length()}, Fisher's Exact test p={pval} ")) +
    scale_y_continuous(breaks = ~ round(unique(pretty(.)))) +
    scale_fill_manual(values = GBMutils::load_color_palette("Verhaak")) +
    guides(fill = guide_legend(title = str_to_title(args$label), nrow = 1, override.aes = list(size = 7))) +
    GBM_theme()
ggsave(p, filename = glue("{args$output_dir}/distribution__{args$label}__{args$method}__{args$suffix}_signatures.pdf"), width = 10, height = 7, dpi = 300)
