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
        description = "TEST-expr-impact",
    )
    parser$add_argument("--input_file", type = "character", help = "Input file")
    parser$add_argument("--min_pct", type = "numeric", help = "Minimum percentage of cells expressing a gene", default = 0.10)
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types/TEST-expr_impact")
    args$meta <- glue("{here::here()}/001_data_local/seurat_annot_adapted__metadata.rds")
    args$min_pct <- 0
    args$input_file <- "output/CellClass_L4_min3_types/100_preprocessing/seurat/6419_cortex.rds"
    args$interactions_db <- glue("{here::here()}/001_data_local/interactions_db/interactions_ref.rds")
    args$interactions <- glue("{here::here()}/output/CellClass_L4_min3_types/400_consensus/401_samples_combined.rds")
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
pacman::p_load(Seurat, ggplot2, ggpubr)

# obj <- readRDS(args$input_file)
interactions <- readRDS(args$interactions)
all_n_interactions_by_sample <- interactions %>%
    select(Sample, source_target, source, target, interaction, Region, is_stringent) %>%
    distinct() %>%
    group_by(is_stringent, Region, Sample, source_target, source, target) %>%
    reframe(n_interactions = n()) %>%
    mutate(setname = case_when(str_detect(source, "Malignant") ~ "Malignant-Other", str_detect(target, "Malignant") ~ "Other-Malignant"))

db <- readRDS(args$interactions_db)
# genes_db <- db %>%
#     select(genename_a, genename_b) %>%
#     pull() %>%
#     unique()
# log_info(glue("Number of genes: {length(genes_db)}"))

# avg_exp <- Seurat::DotPlot(obj, features = rownames(obj), group.by = "CellClass_L4")$data %>% rename(gene = features.plot)
# rownames(avg_exp) <- NULL

# number_of_genes_per_class <- avg_exp %>%
#     filter(pct.exp > args$min_pct, gene %in% genes_db) %>%
#     mutate(dummy = 1) %>%
#     group_by(id) %>%
#     reframe(n = sum(dummy), expressed_genes = paste0(gene, collapse = ", ")) %>%
#     mutate(fraction = n / length(genes_db))

# write.csv(number_of_genes_per_class, glue("{args$output_dir}/{get_name(args$input_file)}__ngenes_per_class.csv"), row.names = FALSE)

# saveRDS(number_of_genes_per_class, glue("{args$output_dir}/{get_name(args$input_file)}__ngenes_per_class.rds"))

meta <- readRDS(args$meta)
rownames(meta) <- NULL

files <- list.files(glue("{args$output_dir}/objects_minpct{args$min_pct*100}"), pattern = "*.rds", full.names = TRUE)

all_obj <- lapply(files, function(file) {
    obj <- readRDS(file)
    name <- str_split(get_name(file), "__", simplify = TRUE)[1]
    obj$name <- name
    return(obj)
})

all_df <- do.call(rbind, all_obj) %>%
    rename(Sample = name) %>%
    left_join(meta %>% select(Sample, Region) %>% distinct(), by = "Sample")
saveRDS(all_df, glue("{args$output_dir}/genes_per_class.rds"))

# log_info(glue("Number of samples: {length(all_df %>% pull(Sample) %>% unique())}"))

# p_fraction <- ggplot(all_df, aes(x = id, y = fraction, color = Region)) +
#     geom_boxplot(outlier.shape = NA) +
#     coord_flip() +
#     custom_theme() +
#     labs(x = "Cell type", title = "Number of ligands/receptors expressed from database per cell type", subtitle = glue("min. fraction of cells that express the ligand/receptor > {args$min_pct}"), y = "Fraction of ligands and receptors expressed") +
#     geom_point(position = position_jitterdodge(), size = .75)

# p_fraction_stat_test <- ggboxplot(all_df, x = "Region", y = "fraction", facet.by = "id", add = "jitter", outlier.shape = NA, ylab = "Fraction of ligands and receptors expressed", xlab = "Cell type", title = "Number of ligands/receptors expressed from database per cell type", subtitle = glue("min. fraction of cells that express the ligand/receptor > {args$min_pct}"), legend = "bottom", ggtheme = theme_pubr(), ncol = 5) + stat_compare_means()
# ggsave(plot = p_fraction_stat_test, filename = glue("{args$output_dir}/fractiongenes_per_class_minpct{args$min_pct * 100}__stat_test.pdf"), width = 10, height = 10, dpi = 300)


# p_fraction
# p_counts <- ggplot(all_df, aes(x = id, y = n, color = Region)) +
#     geom_boxplot(outlier.shape = NA) +
#     coord_flip() +
#     custom_theme() +
#     labs(x = "Cell type", title = "Number of ligands/receptors expressed from database per cell type", subtitle = glue("min. fraction of cells that express the ligand/receptor > {args$min_pct}"), y = "Number of genes") +
#     geom_point(position = position_jitterdodge(), size = .75)

# ggsave(plot = p_fraction, filename = glue("{args$output_dir}/fractiongenes_per_class_minpct{args$min_pct * 100}.pdf"), width = 10, height = 12, dpi = 300)

# ggsave(plot = p_counts, filename =
# glue("{args$output_dir}/ngenes_per_class_minpct{args$min_pct * 100}.pdf"),
# width = 10, height = 12, dpi = 300)

# ---- Number of interactions ---- $
celltypes_oi <- c(all_df %>%
    pull(id) %>%
    unique())
celltype_pairs <- tidyr::crossing(source = celltypes_oi, target = celltypes_oi) %>% mutate(source_target = paste0(source, "__", target))

samples_oi <- all_df %>%
    pull(Sample) %>%
    unique()

n_possible_interactions <- do.call(rbind, lapply(samples_oi, function(sample) {
    df_sub <- all_df %>% filter(Sample == sample)

    celltypes_subset <- celltype_pairs %>% filter(source %in% df_sub$id, target %in% df_sub$id)
    celltypes_subset$Sample <- sample
    n_interactions_by_celltype_pair <- sapply(seq_len(nrow(celltypes_subset)), function(i) {
        source_genes <- df_sub %>%
            filter(id == celltypes_subset[i, "source"] %>% pull()) %>%
            pull(expressed_genes) %>%
            str_split(", ", simplify = TRUE) %>%
            unlist() %>%
            unique()
        target_genes <- df_sub %>%
            filter(id == celltypes_subset[i, "target"] %>% pull()) %>%
            pull(expressed_genes) %>%
            str_split(", ", simplify = TRUE) %>%
            unlist() %>%
            unique()

        n_interactions <- db %>%
            filter(genename_a %in% source_genes, genename_b %in% target_genes) %>%
            pull(interaction) %>%
            unique() %>%
            length()
        return(n_interactions)
    })
    celltypes_subset$n_possible_interactions <- n_interactions_by_celltype_pair
    return(celltypes_subset)
}))

all_n_interactions_by_sample <- all_n_interactions_by_sample %>% left_join(n_possible_interactions %>% select(Sample, source_target, n_possible_interactions))
saveRDS(all_n_interactions_by_sample, glue("{args$output_dir}/all_n_interactions_by_sample.rds"))

for (stringency in c(0, 1)) {
    plt_possible_interactions <- ggplot(data = all_n_interactions_by_sample %>% filter(is_stringent == stringency)) +
        geom_point(aes(x = n_possible_interactions, y = n_interactions, color = target)) +
        facet_grid(Region ~ source) +
        custom_theme() +
        labs(x = "Number of possible interactions", y = "Number of interactions", title = "Number of interactions per cell type pair", subtitle = glue("min. fraction of cells that express the ligand/receptor > {args$min_pct}")) +
        stat_cor(aes(x = n_possible_interactions, y = n_interactions), method = "spearman")

    ggsave(plot = plt_possible_interactions, filename = glue("{args$output_dir}/ninteractions_per_celltypepair_stringency_{stringency}_minpct{args$min_pct * 100}.pdf"), width = 20, height = 8, dpi = 300)


    plot_possible_interactions_overall <- ggplot(data = all_n_interactions_by_sample %>% filter(is_stringent == stringency)) +
        geom_point(aes(x = n_possible_interactions, y = n_interactions)) +
        custom_theme() +
        facet_grid(~Region) +
        labs(x = "Number of possible interactions", y = "Number of interactions", title = "Number of interactions per cell type pair", subtitle = glue("min. fraction of cells that express the ligand/receptor > {args$min_pct}")) +
        stat_cor(aes(x = n_possible_interactions, y = n_interactions), method = "spearman")
    ggsave(plot = plot_possible_interactions_overall, filename = glue("{args$output_dir}/ninteractions_per_celltypepair_stringency_{stringency}_minpct{args$min_pct * 100}__overall.pdf"), width = 10, height = 8, dpi = 300)
}
