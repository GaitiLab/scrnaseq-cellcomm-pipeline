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
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types/TEST-ncells_impact")
    args$percentile <- 0.90
    args$meta <- glue("{here::here()}/001_data_local/seurat_annot_adapted__metadata.rds")
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
pacman::p_load(ggplot2, ggpubr, scales)
meta <- readRDS(args$meta)
interactions <- readRDS(args$interactions)

# Determine number of cells per cell type per sample
n_cells_per_type <- meta %>%
    group_by(Sample, CellClass_L4) %>%
    count()

# Determine number of interactions per pair per sample
all_n_interactions_by_sample <- interactions %>%
    select(Sample, source_target, source, target, interaction, Region, is_stringent) %>%
    distinct() %>%
    group_by(is_stringent, Region, Sample, source_target, source, target) %>%
    reframe(n_interactions = n()) %>%
    mutate(setname = case_when(str_detect(source, "Malignant") ~ "Malignant-Other", str_detect(target, "Malignant") ~ "Other-Malignant")) %>%
    # Add number of cells per cell type for source
    left_join(n_cells_per_type %>% rename(n_cells_source = n), by = c("Sample", "source" = "CellClass_L4")) %>%
    # Add number of cells per cell type for target
    left_join(n_cells_per_type %>% rename(n_cells_target = n), by = c("Sample", "target" = "CellClass_L4"))
saveRDS(all_n_interactions_by_sample, glue("{args$output_dir}/all_n_interactions_by_sample.rds"))

stringency <- 0
for (stringency in c(0, 1)) {
    legend_max <- plyr::round_any(quantile(all_n_interactions_by_sample %>% filter(is_stringent == stringency) %>% pull(n_interactions), args$percentile), 25, f = ceiling)

    model <- lm(n_interactions ~ n_cells_source + n_cells_target, data = all_n_interactions_by_sample %>% filter(is_stringent == stringency))
    r_squared <- sigma(model) / mean(all_n_interactions_by_sample %>% filter(is_stringent == stringency) %>% pull(n_interactions))

    model_sum <- summary(model)
    r_squared <- model_sum$r.squared
    log_info(glue("R^2: {r_squared}"))

    plt <- ggplot(data = all_n_interactions_by_sample %>% filter(is_stringent == stringency)) +
        geom_point(aes(x = n_cells_source, y = n_cells_target, color = n_interactions)) +
        custom_theme() +
        guides(color = guide_colourbar(barwidth = 10, ticks = TRUE)) +
        scale_colour_gradient(
            low = "yellow",
            high = "red",
            space = "Lab",
            na.value = "red",
            guide = "colourbar", aesthetics = "colour", limits = c(1, legend_max), breaks = pretty_breaks()
        ) +
        labs(x = "log10(Number of cells sender)", y = "log10(Number of cells receiver)", title = "Number of interactions", color = "Number of interactions", subtitle = glue("Capped colourbar at {args$percentile * 100}% of interactions overall")) +
        facet_grid(~Region) +
        scale_x_continuous(trans = "log10") +
        scale_y_continuous(trans = "log10")

    ggsave(plot = plt, filename = glue("{args$output_dir}/impact_cell_abundance__stringency_{stringency}.pdf"), width = 15, height = 8)

    # Auto communication
    auto_interactions <- all_n_interactions_by_sample %>% filter(source == target, is_stringent == stringency)
    plt_auto <- ggplot(data = auto_interactions) +
        geom_point(aes(x = n_cells_source, y = n_interactions, color = source, shape = Region)) +
        labs(x = "Number of cells", y = "Number of interactions", title = "Communication of cell group with itself") +
        custom_theme() +
        stat_cor(aes(x = n_cells_source, y = n_interactions), method = "spearman")

    other_combi <- all_n_interactions_by_sample %>%
        filter(source != target, is_stringent == stringency) %>%
        mutate(total_cells = n_cells_source + n_cells_target, product_of_cells = n_cells_source * n_cells_target, max_cells = pmax(n_cells_source, n_cells_target), min_cells = pmin(n_cells_source, n_cells_target), mean_cells = (n_cells_source + n_cells_target) / 2, median_cells = median(c(n_cells_source, n_cells_target)))
    #  Source + target cells
    plt_sum <- ggplot(data = other_combi) +
        geom_point(aes(x = total_cells, y = n_interactions, shape = Region, color = target)) +
        labs(x = "(source + target) cells", y = "Number of interactions", title = "Communication of cell group with other cell group") +
        custom_theme() +
        facet_wrap(~source, nrow = 3) +
        stat_cor(aes(x = total_cells, y = n_interactions), method = "spearman")

    #  Source * target cells
    plt_product <- ggplot(data = other_combi) +
        geom_point(aes(x = product_of_cells, y = n_interactions, shape = Region, color = target)) +
        labs(x = "(source x target) cells", y = "Number of interactions", title = "Communication of cell group with other cell group") +
        custom_theme() +
        facet_wrap(~source, nrow = 3) +
        stat_cor(aes(x = product_of_cells, y = n_interactions), method = "spearman")

    #  Max cells
    plt_max <- ggplot(data = other_combi) +
        geom_point(aes(x = max_cells, y = n_interactions, shape = Region, color = target)) +
        labs(x = "max(source, target) cells", y = "Number of interactions", title = "Communication of cell group with other cell group") +
        custom_theme() +
        facet_wrap(~source, nrow = 3) +
        stat_cor(aes(x = max_cells, y = n_interactions), method = "spearman")

    #  Min cells
    plt_min <- ggplot(data = other_combi) +
        geom_point(aes(x = min_cells, y = n_interactions, shape = Region, color = target)) +
        labs(x = "min(source, target) cells", y = "Number of interactions", title = "Communication of cell group with other cell group") +
        custom_theme() +
        facet_wrap(~source, nrow = 3) +
        stat_cor(aes(x = min_cells, y = n_interactions), method = "spearman")

    #  Mean cells
    plt_mean <- ggplot(data = other_combi) +
        geom_point(aes(x = mean_cells, y = n_interactions, shape = Region, color = target)) +
        labs(x = "mean(source, target) cells", y = "Number of interactions", title = "Communication of cell group with other cell group") +
        custom_theme() +
        facet_wrap(~source, nrow = 3) +
        stat_cor(aes(x = mean_cells, y = n_interactions), method = "spearman")

    all_plts <- list(plt_auto, plt_sum, plt_product, plt_max, plt_min, plt_mean)
    names(all_plts) <- c("plt_auto", "plt_sum", "plt_product", "plt_max", "plt_min", "plt_mean")
    for (i in seq_len(length(all_plts))) {
        print(i)
        ggsave(plot = all_plts[[i]], filename = glue("{args$output_dir}/{names(all_plts)[[i]]}__stringency_{stringency}.pdf"), width = 10, height = 10)
    }
}
