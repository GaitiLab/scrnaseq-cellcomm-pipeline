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
        description = "Inspect relationship between number of cells (sender/receiver) and number of interactions",
    )
    parser$add_argument("--annot", type = "character", help = "Annotation level")
    parser$add_argument("--percentile", type = "numeric", help = "Percentile")
    parser$add_argument("--meta", type = "character", help = "metadata RDS")
    parser$add_argument("--interactions", type = "character", help = "Interactions mvoted (output from 401 script)")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$annot <- "CCI_CellClass_L1"
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/{args$annot}/TEST-ncells_impact")
    args$percentile <- 0.90
    args$meta <- glue("{here::here()}/output/CCI_CellClass_L2/000_data/bw_gbm_regional_study__metadata.rds")
    args$interactions <- glue("{here::here()}/output/{args$annot}/401_combine_samples/401_samples_interactions_mvoted.rds")
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
pacman::p_load(ggplot2, ggpubr, scales, ggtext)

log_info("Load metadata...")
meta <- readRDS(args$meta)

log_info("Load significant interactions per sample...")
interactions <- readRDS(args$interactions) %>% separate(source_target, into = c("source", "target"), sep = "__", remove = FALSE)
# # A tibble: 6 × 16
#   Sample         source_target        source    target    complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb lenient_voting stringent_voting Patient Region_Grouped Batch   Platform
#   <chr>          <chr>                <chr>     <chr>     <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl> <lgl>          <lgl>            <chr>   <chr>          <chr>   <chr>
# 1 6234_2895153_A Malignant__Malignant Malignant Malignant A2M__LRP1                   1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio
# 2 6234_2895153_A Malignant__Malignant Malignant Malignant ACE__BDKRB2                 1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio
# 3 6234_2895153_A Malignant__Malignant Malignant Malignant ADAM10__AXL                 1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio
# 4 6234_2895153_A Malignant__Malignant Malignant Malignant ADAM10__CADM1               1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio
# 5 6234_2895153_A Malignant__Malignant Malignant Malignant ADAM10__GPNMB               1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio
# 6 6234_2895153_A Malignant__Malignant Malignant Malignant ADAM10__IL6R                1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio

# Determine number of cells per cell type per sample
log_info("Determine number of cells per cell type per sample...")
n_cells_per_type <- meta %>%
    group_by(Sample, !!as.symbol(args$annot)) %>%
    count()
# # A tibble: 6 × 3
# # Groups:   Sample, CCI_CellClass_L1 [6]
#   Sample         CCI_CellClass_L1     n
#   <chr>          <chr>            <int>
# 1 6234_2895153_A Astrocyte           19
# 2 6234_2895153_A Endothelial          8
# 3 6234_2895153_A Macrophage           4
# 4 6234_2895153_A Malignant          238
# 5 6234_2895153_A Microglia          103
# 6 6234_2895153_A Neuron             744

log_info("Compute number of interactions per pair per sample...")
types_of_voting <- c("lenient_voting", "stringent_voting")

for (type_of_voting in types_of_voting) {
    log_info(glue("Type of voting: {type_of_voting}..."))
    # Determine number of interactions for each sample x source-target pair + add the number of cells they have for that cell type pair
    all_n_interactions_by_sample <- interactions %>%
        select(Sample, source_target, source, target, complex_interaction, Region_Grouped, !!as.symbol(type_of_voting)) %>%
        distinct() %>%
        filter(!!as.symbol(type_of_voting)) %>%
        group_by(Region_Grouped, Sample, source_target, source, target) %>%
        reframe(n_interactions = n()) %>%
        # Add number of cells per cell type for source
        left_join(n_cells_per_type %>% rename(n_cells_source = n), by = c("Sample", "source" = args$annot)) %>%
        # Add number of cells per cell type for target
        left_join(n_cells_per_type %>% rename(n_cells_target = n), by = c("Sample", "target" = args$annot))
    # # A tibble: 6 × 8
    #   Region_Grouped Sample         source_target              source    target          n_interactions n_cells_source n_cells_target
    #   <chr>          <chr>          <chr>                      <chr>     <chr>                    <int>          <int>          <int>
    # 1 PT             6234_2895153_B Malignant__Malignant       Malignant Malignant                  156            288            288
    # 2 PT             6234_2895153_B Malignant__Microglia       Malignant Microglia                  174            288            100
    # 3 PT             6234_2895153_B Malignant__Neuron          Malignant Neuron                      89            288            540
    # 4 PT             6234_2895153_B Malignant__Oligodendrocyte Malignant Oligodendrocyte             97            288            245
    # 5 PT             6234_2895153_B Microglia__Malignant       Microglia Malignant                  143            100            288
    # 6 PT             6234_2895153_B Microglia__Microglia       Microglia Microglia                  258            100            100

    saveRDS(all_n_interactions_by_sample, glue("{args$output_dir}/all_n_interactions_by_sample_{type_of_voting}.rds"))

    legend_max <- plyr::round_any(quantile(all_n_interactions_by_sample %>% pull(n_interactions), args$percentile), 25, f = ceiling)

    model <- lm(n_interactions ~ n_cells_source + n_cells_target, data = all_n_interactions_by_sample)
    r_squared <- sigma(model) / mean(all_n_interactions_by_sample %>% pull(n_interactions))

    model_sum <- summary(model)
    r_squared <- model_sum$r.squared
    log_info(glue("R^2: {r_squared}"))

    plt <- ggplot(data = all_n_interactions_by_sample) +
        geom_point(aes(x = n_cells_source, y = n_cells_target, color = n_interactions)) +
        custom_theme() +
        guides(color = guide_colourbar(barwidth = 10, ticks = TRUE)) +
        scale_colour_gradient(
            low = "yellow",
            high = "red",
            space = "Lab",
            # Capped at args$percentile (90%), so these values would be considerd NA, therefore manually setting to max (red)
            na.value = "red",
            guide = "colourbar", aesthetics = "colour", limits = c(1, legend_max), breaks = pretty_breaks()
        ) +
        labs(
            x = "Number of cells sender (log<sub>10</sub>)", y = "Number of cells receiver (log<sub>10</sub>)",
            title = "Relationship between number of cells and detected number of interactions",
            color = "Number of interactions",
            subtitle = glue("Capped colourbar at {args$percentile * 100}% of interactions overall")
        ) +
        facet_grid(~Region_Grouped) +
        # scale_x_continuous(trans = "log10") +
        # scale_y_continuous(trans = "log10")
        scale_x_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
    plt
    ggsave(plot = plt, filename = glue("{args$output_dir}/impact_cell_abundance__{type_of_voting}.pdf"), width = 15, height = 8)

    # # Auto communication
    # auto_interactions <- all_n_interactions_by_sample %>% filter(source == target)
    # # plt_auto <- ggplot(data = auto_interactions) +
    # #     geom_point(aes(x = n_cells_source, y = n_interactions, color = source, shape = Region_Grouped)) +
    # #     labs(x = "Number of cells", y = "Number of interactions", title = "Communication of cell group with itself") +
    # #     custom_theme() +
    # #     stat_cor(aes(x = n_cells_source, y = n_interactions), method = "spearman") +  scale_x_log10(
    # #          breaks = scales::trans_breaks("log10", function(x) 10^x),
    # #          labels = scales::trans_format("log10", scales::math_format(10^.x))
    # #      )

    # ggscatter(data = auto_interactions, x = "n_cells_source", y = "n_interactions", color = "source", add = "reg.line", conf.int = TRUE) +
    # stat_cor(method = "spearman", aes(color = source), label.x = 10)

    # # + scale_x_log10(
    # #     breaks = scales::trans_breaks("log10", function(x) 10^x),
    # #     labels = scales::trans_format("log10", scales::math_format(10^.x))
    # # )

    # other_combi <- all_n_interactions_by_sample %>%
    #     filter(source != target, is_stringent == stringency) %>%
    #     mutate(total_cells = n_cells_source + n_cells_target, product_of_cells = n_cells_source * n_cells_target, max_cells = pmax(n_cells_source, n_cells_target), min_cells = pmin(n_cells_source, n_cells_target), mean_cells = (n_cells_source + n_cells_target) / 2, median_cells = median(c(n_cells_source, n_cells_target)))
    # #  Source + target cells
    # plt_sum <- ggplot(data = other_combi) +
    #     geom_point(aes(x = total_cells, y = n_interactions, shape = Region_Grouped, color = target)) +
    #     labs(x = "(source + target) cells", y = "Number of interactions", title = "Communication of cell group with other cell group") +
    #     custom_theme() +
    #     facet_wrap(~source, nrow = 3) +
    #     stat_cor(aes(x = total_cells, y = n_interactions), method = "spearman")

    # #  Source * target cells
    # plt_product <- ggplot(data = other_combi) +
    #     geom_point(aes(x = product_of_cells, y = n_interactions, shape = Region_Grouped, color = target)) +
    #     labs(x = "(source x target) cells", y = "Number of interactions", title = "Communication of cell group with other cell group") +
    #     custom_theme() +
    #     facet_wrap(~source, nrow = 3) +
    #     stat_cor(aes(x = product_of_cells, y = n_interactions), method = "spearman")

    # #  Max cells
    # plt_max <- ggplot(data = other_combi) +
    #     geom_point(aes(x = max_cells, y = n_interactions, shape = Region_Grouped, color = target)) +
    #     labs(x = "max(source, target) cells", y = "Number of interactions", title = "Communication of cell group with other cell group") +
    #     custom_theme() +
    #     facet_wrap(~source, nrow = 3) +
    #     stat_cor(aes(x = max_cells, y = n_interactions), method = "spearman")

    # #  Min cells
    # plt_min <- ggplot(data = other_combi) +
    #     geom_point(aes(x = min_cells, y = n_interactions, shape = Region_Grouped, color = target)) +
    #     labs(x = "min(source, target) cells", y = "Number of interactions", title = "Communication of cell group with other cell group") +
    #     custom_theme() +
    #     facet_wrap(~source, nrow = 3) +
    #     stat_cor(aes(x = min_cells, y = n_interactions), method = "spearman")

    # #  Mean cells
    # plt_mean <- ggplot(data = other_combi) +
    #     geom_point(aes(x = mean_cells, y = n_interactions, shape = Region_Grouped, color = target)) +
    #     labs(x = "mean(source, target) cells", y = "Number of interactions", title = "Communication of cell group with other cell group") +
    #     custom_theme() +
    #     facet_wrap(~source, nrow = 3) +
    #     stat_cor(aes(x = mean_cells, y = n_interactions), method = "spearman")

    # all_plts <- list(plt_auto, plt_sum, plt_product, plt_max, plt_min, plt_mean)
    # names(all_plts) <- c("plt_auto", "plt_sum", "plt_product", "plt_max", "plt_min", "plt_mean")
    # for (i in seq_len(length(all_plts))) {
    #     print(i)
    #     ggsave(plot = all_plts[[i]], filename = glue("{args$output_dir}/{names(all_plts)[[i]]}__stringency_{stringency}.pdf"), width = 10, height = 10)
    # }
}
