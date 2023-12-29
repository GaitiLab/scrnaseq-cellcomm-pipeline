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
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types/503b_number_of_interactions_sample_wise")
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

# ---- Boxplots w/ number of interactions per sample ---- #
n_interactions_by_sample <- interactions %>%
    select(Sample, source_target, source, target, interaction, Region, is_stringent) %>%
    distinct() %>%
    group_by(is_stringent, Region, Sample, source_target, source, target) %>%
    reframe(n_interactions = n()) %>%
    mutate(setname = case_when(str_detect(source, "Malignant") ~ "Malignant-Other", str_detect(target, "Malignant") ~ "Other-Malignant")) %>%
    mutate(Region = factor(Region, levels = c("PT", "TE", "SC")))

for (stringency in c(0, 1)) {
    # Malignant (sender/ligand) to other cell type (receiver/receptor)
    boxplot1 <- n_interactions_by_sample %>%
        filter(
            setname == "Malignant-Other",
            is_stringent == stringency
        ) %>%
        ggplot(data = .) +
        custom_theme() +
        geom_boxplot(aes(x = reorder(target, n_interactions, mean, order = TRUE), y = n_interactions, fill = Region), outlier.shape = NA) +
        facet_grid(~Region, scales = "free_x") +
        geom_jitter(aes(x = reorder(target, n_interactions, mean, order = TRUE), y = n_interactions), size = 0.4) +
        rotate_x_text() +
        labs(x = "Receiver", y = "Number of interactions", title = "Malignant (sender w/ ligand) to other cell type (receiver w/ receptor)")

    # Other cell type (sender/ligand) to Malignant (receiver/receptor)
    boxplot2 <- n_interactions_by_sample %>%
        filter(
            setname == "Other-Malignant",
            is_stringent == stringency
        ) %>%
        ggplot(data = .) +
        geom_boxplot(aes(x = reorder(source, n_interactions, mean, order = TRUE), y = n_interactions, fill = Region), outlier.shape = NA) +
        facet_grid(~Region, scales = "free_x") +
        geom_jitter(aes(x = reorder(source, n_interactions, mean, order = TRUE), y = n_interactions), size = 0.4) +
        custom_theme() +
        rotate_x_text() +
        labs(x = "Sender", y = "Number of interactions", title = "Other cell type (sender w/ ligand) to Malignant (receiver w/ receptor)")

    box_combined <- ggarrange(boxplot1, boxplot2, nrow = 2)
    ggsave(plot = box_combined, filename = glue("{args$output_dir}/boxplot__stringency_{stringency}.pdf"), width = 10, height = 10)
}

df <- interactions %>%
    select(is_stringent, Sample, Region, source_target, interaction, source, target) %>%
    group_by(is_stringent, Sample, Region, source_target, source, target) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    select(is_stringent, Sample, Region, source, target, n) %>%
    pivot_wider(names_from = target, values_from = n, values_fill = 0)

df_norm <- df %>% mutate(across(c(interactions %>% pull(target) %>% unique()), function(x) {
    x / n_possible_interactions
}))



# sample <- available_samples[1]
# stringency <- 0

# Create heatmaps for each individual sample
for (sample in available_samples) {
    for (stringency in c(0, 1)) {
        # Determine max value for heatmap legends
        max_fraction <- df_norm %>%
            filter(is_stringent == stringency) %>%
            select(all_of(c(interactions %>% pull(target) %>% unique()))) %>%
            max()
        legend_norm_max <- plyr::round_any(max_fraction, 0.05, f = ceiling)

        max_n <- df %>%
            filter(is_stringent == stringency) %>%
            select(all_of(c(interactions %>% pull(target) %>% unique()))) %>%
            max()
        legend_max <- plyr::round_any(max_n, 50, f = ceiling)

        region <- df %>%
            filter(Sample == sample) %>%
            pull(Region) %>%
            unique()
        log_info(glue("Region: {region}"))

        mat <- df %>%
            filter(is_stringent == stringency, Sample == sample, Region == region) %>%
            ungroup() %>%
            select(-Region, -is_stringent, -Sample) %>%
            ungroup() %>%
            column_to_rownames("source") %>%
            data.matrix()

        mat <- mat[order(rownames(mat)), order(colnames(mat))]

        create_hm(mat = mat, col_fun = circlize::colorRamp2(c(0, legend_max), c("white", "red")), output_file = glue("{args$output_dir}/heatmap_{sample}__{region}__stringency_{stringency}.pdf"), legend_title = "Interactions", save_plot = TRUE, custom_cell_fun = custom_cell_function, column_title = "Receiver", row_title = "Sender", column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE)

        mat_norm <- df_norm %>%
            filter(is_stringent == stringency, Sample == sample, Region == region) %>%
            ungroup() %>%
            select(-Region, -is_stringent, -Sample) %>%
            ungroup() %>%
            column_to_rownames("source") %>%
            data.matrix()

        mat_norm <- mat_norm[order(rownames(mat_norm)), order(colnames(mat_norm))]

        create_hm(mat = mat_norm, col_fun = circlize::colorRamp2(c(0, legend_norm_max), c("white", "red")), output_file = glue("{args$output_dir}/heatmap_{sample}__{region}__stringency_{stringency}__norm.pdf"), legend_title = "Fraction detected", save_plot = TRUE, custom_cell_fun = custom_cell_function_default, column_title = "Receiver", row_title = "Sender", column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE)
    }
}

log_info("COMPLETED!")
