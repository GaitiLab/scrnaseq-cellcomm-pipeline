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
        description = "Comparison runs",
    )

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/comparing_runs")
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
run_name1 <- "CCI_CellClass_L2_conf_malign"
run_name2 <- "CCI_CellClass_L2_conf_min50"

# Setup grid
masked <- c("lenient_condition", "stringent_condition")
source_targets_oi <- readRDS(glue("{here::here()}/output/{run_name1}/402_aggregation/402_patient_interactions_mvoted_w_filters.rds")) %>%
    pull(source_target) %>%
    unique()

combi_grid <- data.frame(expand.grid(masked, source_targets_oi, REGION_GROUPED_LEVELS))
colnames(combi_grid) <- c("mask", "source_target", "region")

comparison <- pbapply::pblapply(seq_len(nrow(combi_grid)), function(row_id) {
    mask <- as.character(combi_grid[row_id, "mask"])
    current_source_target <- as.character(combi_grid[row_id, "source_target"])
    current_region <- as.character(combi_grid[row_id, "region"])

    group1 <- readRDS(glue("{here::here()}/output/{run_name1}/402_aggregation/402_patient_interactions_mvoted_w_filters.rds")) %>% filter(!!as.symbol(mask))
    group2 <- readRDS(glue("{here::here()}/output/{run_name2}/402_aggregation/402_patient_interactions_mvoted_w_filters.rds")) %>% filter(!!as.symbol(mask))

    group1 <- group1 %>%
        filter(source_target == current_source_target, Region_Grouped == current_region) %>%
        pull(complex_interaction)

    group2 <- group2 %>%
        filter(source_target == current_source_target, Region_Grouped == current_region) %>%
        pull(complex_interaction)

    common <- intersect(group1, group2)
    overlap_frac_wrt_union <- length(common) / length(union(group1, group2))
    overlap_frac_wrt_group1 <- length(common) / length(group1)
    overlap_frac_wrt_group2 <- length(common) / length(group2)

    return(c(
        source_target = current_source_target,
        region = current_region,
        filter_option = mask,
        group1 = length(group1),
        group2 = length(group2),
        n_common = length(common),
        overlap_frac_wrt_union = overlap_frac_wrt_union,
        overlap_frac_wrt_group2 = overlap_frac_wrt_group2,
        overlap_frac_wrt_group1 = overlap_frac_wrt_group1
    ))
})

comparison_df <- data.frame(do.call(rbind, comparison)) %>%
    mutate(
        group1 = as.numeric(group1),
        group2 = as.numeric(group2),
        n_common = as.numeric(n_common),
        overlap_frac_wrt_union = as.numeric(overlap_frac_wrt_union),
        overlap_frac_wrt_group1 = as.numeric(overlap_frac_wrt_group1),
        overlap_frac_wrt_group2 = as.numeric(overlap_frac_wrt_group2),
        filter_option = str_replace_all(filter_option, "_condition", "")
    ) %>%
    filter(group1 > 0, group2 > 0)
colnames(comparison_df)[colnames(comparison_df) %in% c("overlap_frac_wrt_group1", "overlap_frac_wrt_group2")] <- c(paste0("overlap_wrt\n", run_name1), paste0("overlap_wrt\n", run_name2))


write.csv(comparison_df, file = glue("{args$output_dir}/comp_{run_name1}_vs_{run_name2}.csv"))

comparison_df <- comparison_df %>%
    melt(id = c("group1", "group2", "n_common", "source_target", "region", "filter_option"))
head(comparison_df)
# smallest_overlap

overlap_options <- comparison_df %>%
    pull(variable) %>%
    unique()

# smallest_overlap <- comparison_df %>%
#     pull(value) %>%
#     min()

log_info("Visualize...")
hist_plt <- ggplot(data = comparison_df) +
    geom_histogram(aes(x = value)) +
    ggh4x::facet_nested(region + filter_option ~ variable) +
    custom_theme() +
    labs(
        x = "Fraction of overlap", y = "Counts", title = glue("Comparison: {run_name1} vs {run_name2}"),
        # subtitle = glue("Smallest frac of overlap: {smallest_overlap}")
    ) +
    geom_vline(xintercept = 0.7, color = "red", linetype = "dashed")
ggsave(plot = hist_plt, filename = glue("{args$output_dir}/comp_{run_name1}_vs_{run_name2}.pdf"), width = 15, height = 15)
