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
        description = "Visualize number of interactions per sample",
    )
    parser$add_argument("--input_file", type = "character", help = "Path to input file")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types/503a_number_of_interactions_union")
    args$input_file <- glue("{here::here()}/output/CellClass_L4_min3_types/400_consensus/401_samples_combined.rds")
    args$interactions_db <- "001_data_local/interactions_db/interactions_ref.rds"
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
pacman::p_load(ggplot2, ggpubr)
interactions <- readRDS(args$input_file)

available_samples <- interactions %>%
    select(Sample) %>%
    distinct() %>%
    pull(Sample) %>%
    unique()
log_info(glue("Number of samples: {length(available_samples)}"))


# Load interactions database
log_info("Determine number of possible interactions...")
interactions_db <- readRDS(args$interactions_db)
n_possible_interactions <- nrow(interactions_db)
log_info(glue("Number of possible interactions: {n_possible_interactions}"))

# Combine all samples (no-filtering -> taking the union)
df <- interactions %>%
    # Check for distinct interactions regardless of sample, identified in at
    # least 1 sample for a region
    select(is_stringent, Region, source_target, interaction, source, target) %>%
    distinct() %>%
    group_by(is_stringent, Region, source_target, source, target) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    select(is_stringent, Region, source, target, n) %>%
    pivot_wider(names_from = target, values_from = n, values_fill = 0)

df_norm <- df %>% mutate(across(c(interactions %>% pull(target) %>% unique()), function(x) {
    x / n_possible_interactions
}))


# Number of interactions by pair
n_interactions_by_pair <- interactions %>%
    select(is_stringent, Region, source_target, source, target, interaction) %>%
    distinct() %>%
    group_by(is_stringent, Region, source_target, source, target) %>%
    summarise(n_total = n())

# Check how many interactions are detected in > 1 samples.
n_samples_by_interactions <- interactions %>%
    select(is_stringent, Sample, Region, source_target, source, target, interaction) %>%
    distinct() %>%
    group_by(is_stringent, Region, source_target, source, target, interaction) %>%
    summarise(n_samples = n())

# How many samples have > 1 interaction for a given interaction pair?
n_samples_by_interactions <- n_samples_by_interactions %>%
    mutate(multi_samples = n_samples > 1) %>%
    group_by(is_stringent, Region, source_target, source, target) %>%
    summarise(n_detected = sum(multi_samples))

combined <- merge(n_interactions_by_pair, n_samples_by_interactions, by = c("is_stringent", "Region", "source_target", "source", "target")) %>% mutate(fraction = n_detected / n_total, n_support = paste0(n_detected, "/", n_total))

# ggplot() +
#     geom_histogram(data = combined, aes(x = fraction), bins = 100) +
#     facet_grid(is_stringent ~ Region) +
#     custom_theme()

# ggplot() +
#     geom_violin(data = combined, aes(x = Region, y = fraction, fill = is_stringent)) +
#     custom_theme()

available_regions <- df %>%
    pull(Region) %>%
    unique()

for (stringency in c(0, 1)) {
    # Determine max value for heatmap legends
    max_fraction <- df_norm %>%
        filter(is_stringent == stringency) %>%
        select(all_of(c(interactions %>% pull(target) %>% unique()))) %>%
        max()
    legend_norm_max <- plyr::round_any(max_fraction, 0.025, f = ceiling)

    max_n <- df %>%
        filter(is_stringent == stringency) %>%
        select(all_of(c(interactions %>% pull(target) %>% unique()))) %>%
        max()
    legend_max <- plyr::round_any(max_n, 50, f = ceiling)

    for (region in available_regions) {
        log_info(glue("\n\nVisualizing: \nstringency={stringency}, region={region}..."))

        combi_df <- combined %>% filter(is_stringent == stringency, Region == region)
        plt <- ggplot(data = combi_df, aes(x = source, y = target)) +
            geom_tile(aes(fill = fraction), colour = "white", size = 2) +
            custom_theme() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            labs(x = "Sender", y = "Region", fill = "Fraction identified\n> 1 sample") +
            guides(fill = guide_colourbar(barwidth = 10, label.position = "bottom")) +
            custom_theme() +
            geom_text(aes(label = n_support)) +
            scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5, limits = c(0, 1))
        ggsave(plot = plt, filename = glue("{args$output_dir}/heatmap__hard_combine__{region}__stringency_{stringency}_support.pdf"), width = 10, height = 10)

        mat <- df %>%
            filter(Region == region, is_stringent == stringency) %>%
            select(-Region, -is_stringent) %>%
            ungroup() %>%
            column_to_rownames("source") %>%
            data.matrix()


        # remove rows and columns that are all 0
        mat <- mat[, which(colSums(mat) != 0)]
        mat <- mat[which(rowSums(mat) != 0), ]

        mat <- mat[order(rownames(mat)), order(colnames(mat))]

        create_hm(mat = mat, col_fun = circlize::colorRamp2(c(0, legend_max), c("white", "red")), output_file = glue("{args$output_dir}/heatmap__hard_combine__{region}__stringency_{stringency}.pdf"), legend_title = "Interactions", save_plot = TRUE, custom_cell_fun = custom_cell_function, column_title = "Receiver", row_title = "Sender", column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE)

        mat_norm <- df_norm %>%
            filter(Region == region, is_stringent == stringency) %>%
            select(-Region, -is_stringent) %>%
            ungroup() %>%
            column_to_rownames("source") %>%
            data.matrix()

        # remove rows and columns that are all 0
        mat_norm <- mat_norm[, which(colSums(mat_norm) != 0)]
        mat_norm <- mat_norm[which(rowSums(mat_norm) != 0), ]

        mat_norm <- mat_norm[order(rownames(mat_norm)), order(colnames(mat_norm))]


        create_hm(mat = mat_norm, col_fun = circlize::colorRamp2(c(0, legend_norm_max), c("white", "red")), output_file = glue("{args$output_dir}/heatmap__hard_combine__{region}__stringency_{stringency}__norm.pdf"), legend_title = "Interactions", save_plot = TRUE, custom_cell_fun = custom_cell_function_default, column_title = "Receiver", row_title = "Sender", column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE)
    }
}

log_info("COMPLETED!")
