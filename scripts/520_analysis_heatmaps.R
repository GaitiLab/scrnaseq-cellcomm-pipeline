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
        description = "Visualize presence of interactions",
    )
    parser$add_argument("--all_interactions", type = "character", help = "Path to all interactions")
    parser$add_argument("--interactions_oi", type = "character", help = "Path to interactions of interest")
    parser$add_argument("--is_stringent", type = "integer", help = "Whether to use stringent (all methods) or relaxed interactions (LIANA + other method)", default = 1)
    parser$add_argument("--literature_findings", type = "character", default = glue("{here::here()}/data/literature_findings_v3.rds"))
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types/500_analysis_heatmaps")
    args$interactions_oi <- glue("{here::here()}/output/CellClass_L4_min3_types/400_consensus/400c_post_filtering.rds")
    args$all_interactions <- glue("{here::here()}/output/CellClass_L4_min3_types/400_consensus/401_samples_combined.rds")
    args$literature_findings <- glue("{here::here()}/001_data_local/literature_findings_v3.rds")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- paste0(args$output_dir)
create_dir(output_dir)

# Load additional libraries
pacman::p_load_gh("jokergoo/ComplexHeatmap", "jmw86069/colorjam")
pacman::p_load(circlize)

log_info(glue("\n\nParameters: \nInteraction file:\t{args$all_interactions}\nCCIs of interest:\t{args$interactions_oi}\nStringency:\t\t{args$is_stringent}"))

log_info("Load interactions...")
all_interactions <- readRDS(args$all_interactions) %>%
    group_by(is_stringent, interaction, Region_Grouped, source_target, source, target) %>%
    reframe(
        n_samples = n(),
        samples = paste0(Sample, collapse = ", "),
        pval_combined = get_p(pval)
    ) %>%
    mutate(setname = case_when(
        source == "Malignant" & target == "Malignant" ~ "Malignant-Malignant",
        str_detect(source, "Malignant") ~ "Malignant-Other", str_detect(target, "Malignant") ~ "Other-Malignant"
    )) %>%
    mutate(Region_Grouped = factor(Region_Grouped, levels = c("PT", "TE", "SC", "NC")))

log_info("Load CCIs of interest (post-filtered)...")
interactions_oi <- readRDS(args$interactions_oi)

available_regions <- interactions_oi %>%
    pull(Region_Grouped) %>%
    unique()

# region <- available_regions[1]
# stringency <- 1

for (stringency in c(0, 1)) {
    for (region in available_regions) {
        log_info(glue("Stringency: {stringency}, Region_Grouped: {region}"))
        # Malignant-Other
        setname_oi <- "Malignant-Other"
        interactions_oi_subset <- interactions_oi %>%
            filter(setname == setname_oi, Region_Grouped == region, cond_min_samples_region, is_stringent == stringency) %>%
            pull(interaction) %>%
            unique()

        mat <- all_interactions %>%
            filter(interaction %in% interactions_oi_subset, setname == setname_oi, Region_Grouped == region, is_stringent == stringency) %>%
            select(target, interaction, n_samples) %>%
            pivot_wider(names_from = target, values_from = n_samples, values_fill = 0) %>%
            column_to_rownames("interaction") %>%
            as.matrix()

        mat_pval <- all_interactions %>%
            filter(interaction %in% interactions_oi_subset, setname == setname_oi, Region_Grouped == region, is_stringent == stringency) %>%
            select(target, interaction, pval_combined) %>%
            pivot_wider(names_from = target, values_from = pval_combined, values_fill = NA) %>%
            column_to_rownames("interaction") %>%
            as.matrix()

        is_detected <- rowSums(mat > 0)
        total_samples <- rowSums(mat)

        row_order <- order(is_detected, decreasing = TRUE)
        mat_sorted <- mat[row_order, ]
        mat_pval_sorted <- mat_pval[row_order, ]

        is_detected <- colSums(mat_sorted > 0)
        total_samples <- colSums(mat_sorted)

        col_order <- order(is_detected, decreasing = TRUE)
        mat_sorted <- mat_sorted[, col_order]
        mat_pval_sorted <- mat_pval_sorted[, col_order]

        max_vote <- plyr::round_any(max(mat_sorted), 2, f = ceiling)
        col_fun <- colorRamp2(c(0, 1, max_vote), c("white", "yellow", "red"))

        cell_fun <- function(j, i, x, y, width, height, fill) {
            grid.rect(
                x = x, y = y, width = width, height = height,
                gp = gpar(col = "grey", fill = "white")
            )
            if (!is.na(mat_pval_sorted[i, j])) {
                grid.circle(x = x, y = y, r = -log10(mat_pval_sorted[i, j]) / 10 * min(unit.c(width, height)), gp = gpar(fill = col_fun(mat_sorted[i, j]), col = NA))
            }
        }
        create_hm(mat = mat_sorted, output_file = glue("{output_dir}/heatmap__{setname_oi}__{region}__stringency_{stringency}__pval.pdf"), legend_title = "Samples", save_plot = TRUE, column_title = "Receiver", row_title = "Interaction", column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE, custom_cell_fun = cell_fun, col_fun = col_fun, cell_size = 4, heatmap_legend_list = list(
            Legend(
                labels = c("lower", " ", " ", " ", "higher"),
                title = "\n-log10(p-value)",
                type = "points",
                # ncol = 5,
                pch = 16,
                size = unit(1:5, "mm"),
                direction = "horizontal",
                legend_gp = gpar(
                    col = "black"
                ),
                background = "white"
            )
        ))

        create_hm(mat = mat_sorted, output_file = glue("{output_dir}/heatmap__{setname_oi}__{region}__stringency_{stringency}.pdf"), legend_title = "Samples", save_plot = TRUE, column_title = "Receiver", row_title = "Interaction", column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE, custom_layer_fun = custom_layer_function, col_fun = col_fun, cell_size = 4, )


        # Other-Malignant
        setname_oi <- "Other-Malignant"
        interactions_oi_subset <- interactions_oi %>%
            filter(setname == setname_oi, Region_Grouped == region, cond_min_samples_region, is_stringent == stringency) %>%
            pull(interaction) %>%
            unique()

        mat <- all_interactions %>%
            filter(interaction %in% interactions_oi_subset, setname == setname_oi, Region_Grouped == region, is_stringent == stringency) %>%
            select(source, interaction, n_samples) %>%
            pivot_wider(names_from = source, values_from = n_samples, values_fill = 0) %>%
            column_to_rownames("interaction") %>%
            as.matrix()

        mat_pval <- all_interactions %>%
            filter(interaction %in% interactions_oi_subset, setname == setname_oi, Region_Grouped == region, is_stringent == stringency) %>%
            select(source, interaction, pval_combined) %>%
            pivot_wider(names_from = source, values_from = pval_combined, values_fill = NA) %>%
            column_to_rownames("interaction") %>%
            as.matrix()

        is_detected <- rowSums(mat > 0)
        total_samples <- rowSums(mat)

        row_order <- order(is_detected, decreasing = TRUE)
        mat_sorted <- mat[row_order, ]
        mat_pval_sorted <- mat_pval[row_order, ]

        is_detected <- colSums(mat_sorted > 0)
        total_samples <- colSums(mat_sorted)

        col_order <- order(is_detected, decreasing = TRUE)
        mat_sorted <- mat_sorted[, col_order]
        mat_pval_sorted <- mat_pval_sorted[, col_order]

        max_vote <- plyr::round_any(max(mat_sorted), 2, f = ceiling)
        col_fun <- colorRamp2(c(0, 1, max_vote), c("white", "yellow", "red"))

        cell_fun <- function(j, i, x, y, width, height, fill) {
            grid.rect(
                x = x, y = y, width = width, height = height,
                gp = gpar(col = "grey", fill = "white")
            )
            if (!is.na(mat_pval_sorted[i, j])) {
                grid.circle(x = x, y = y, r = -log10(mat_pval_sorted[i, j]) / 10 * min(unit.c(width, height)), gp = gpar(fill = col_fun(mat_sorted[i, j]), col = NA))
            }
        }
        create_hm(mat = mat_sorted, output_file = glue("{output_dir}/heatmap__{setname_oi}__{region}__stringency_{stringency}__pval.pdf"), legend_title = "Samples", save_plot = TRUE, column_title = "Sender", row_title = "Interaction", column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE, custom_cell_fun = cell_fun, col_fun = col_fun, cell_size = 4, heatmap_legend_list = list(
            Legend(
                labels = c("lower", " ", " ", " ", "higher"),
                title = "\n-log10(p-value)",
                type = "points",
                # ncol = 5,
                pch = 16,
                size = unit(1:5, "mm"),
                direction = "horizontal",
                legend_gp = gpar(
                    col = "black"
                ),
                background = "white"
            )
        ))

        create_hm(mat = mat_sorted, output_file = glue("{output_dir}/heatmap__{setname_oi}__{region}__stringency_{stringency}.pdf"), legend_title = "Samples", save_plot = TRUE, column_title = "Sender", row_title = "Interaction", column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE, custom_layer_fun = custom_layer_function, col_fun = col_fun, cell_size = 4, )
    }
}
