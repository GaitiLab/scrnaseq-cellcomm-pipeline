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
        description = "Simple heatmaps representing number of interactions between cell type groups"
    )
    parser$add_argument("--condition_varname", type = "character", default = "Region_Grouped")
    parser$add_argument("--annot", type = "character", default = "CCI_CellClass_L1")
    parser$add_argument("--input_file", type = "character", help = "Aggregated interactions dataframe (RDS file)")
    parser$add_argument("--agg_level", type = "character", help = "Level of aggregation: sample or patient")

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    run_name <- "CCI_CellClass_L1_w_agg"
    args$annot <- "CCI_CellClass_L1"
    args$output_dir <- glue("{here::here()}/output/{run_name}/501_heatmap_n_interactions")
    args$input_file <- glue("{here::here()}/output/{run_name}/402_aggregation/402_interactions_combi_agg.rds")
    args$condition_varname <- "Region_Grouped"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- glue("{args$output_dir}/")
create_dir(output_dir)

# Load additional libraries

log_info("Load interactions...")
interactions <- readRDS(args$input_file) %>% filter(pval < 0.05)

# Count the number of interactions between cell type groups
# Types of filters:
# - stringent_region (voting method based stringent + take into account only region)
# - stringent_region_pair (voting stringent + take into account both region and presence of pair in sample of that region)
# - lenient_region (take into account only region)
# - lenient_region_pair (take into account both region and presence of pair in sample of that region)
# INTERACTIONS_POST_FILTERING_OPTIONS <- c("stringent_condition", "stringent_condition_pair", "lenient_condition", "lenient_condition_pair")
avail_regions <- interactions %>%
    pull(!!sym(args$condition_varname)) %>%
    unique()

for (option in INTERACTIONS_POST_FILTERING_OPTIONS) {
    if (!option %in% colnames(interactions)) {
        next
    }
    for (current_region in avail_regions) {
        # option <- "lenient_region"
        # current_region <- "PT"

        interactions_filtered <- interactions %>%
            ungroup() %>%
            filter(!!sym(option), !is.na(!!sym(option))) %>%
            select(!!sym(args$condition_varname), source_target, complex_interaction) %>%
            group_by(!!sym(args$condition_varname), source_target) %>%
            separate(source_target, c("source", "target"), sep = "__") %>%
            group_by(!!sym(args$condition_varname), source, target) %>%
            summarise(n = n()) %>%
            ungroup()

        log_info("Current selection...")
        interactions_subset_long <- interactions_filtered %>%
            filter(!!sym(args$condition_varname) == current_region) %>%
            ungroup() %>%
            select(source, target, n) %>%
            # mutate(n = n**2) %>%
            mutate(n_norm = (n - min(n)) / (max(n) - min(n)))

        log_info("Create heatmap w/ matrix...")
        interactions_wide <- interactions_subset_long %>%
            select(source, target, n) %>%
            pivot_wider(names_from = target, values_from = n, values_fill = 0) %>%
            column_to_rownames("source")

        mat <- data.matrix(interactions_wide)
        mat <- mat[, which(colSums(mat) != 0)]
        mat <- mat[which(rowSums(mat) != 0), ]
        mat <- mat[order(rownames(mat)), order(colnames(mat))]
        legend_max <- plyr::round_any(max(mat), 25, f = ceiling)

        create_hm(
            mat = mat, col_fun = circlize::colorRamp2(c(0, legend_max), c("white", "red")),
            output_file = glue("{output_dir}/heatmap__{option}__{current_region}.pdf"),
            legend_title = "Interactions", save_plot = TRUE,
            custom_cell_fun = custom_cell_function,
            column_title = "Receiver", row_title = "Sender",
            column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE
        )
        # log_info(glue("Option: {option}, Region: {current_region}"))
        # log_info("Malignant as receiver:")
        # # print(sort(mat[, "Malignant"], decreasing = TRUE))
        # # log_info("Malignant as sender:")
        # # print(sort(mat["Malignant", ], decreasing = TRUE))

        # log_info("Create heatmap w/ malignant axes only...")
        # malignant_target <- mat[, "Malignant"]
        # malignant_sender <- mat["Malignant", ]

        # mat <- rbind(malignant_sender, malignant_target)
        # rownames(mat) <- c("Sender", "Receiver")

        # legend_norm_max <- plyr::round_any(max(mat), 25, f = ceiling)
        # create_hm(
        #     mat = mat, col_fun = circlize::colorRamp2(c(0, legend_norm_max), c("white", "red")),
        #     output_file = glue("{args$output_dir}/heatmap_malignant_axis_only__{option}__{current_region}.pdf"),
        #     legend_title = "Interactions", save_plot = TRUE,
        #     custom_cell_fun = custom_cell_function,
        #     column_title = "", row_title = "",
        #     column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE
        # )


        # log_info("Create heatmap w/ normalized matrix...")
        # interactions_norm_wide <- interactions_subset_long %>%
        #     select(source, target, n_norm) %>%
        #     pivot_wider(names_from = target, values_from = n_norm, values_fill = 0) %>%
        #     column_to_rownames("source")

        # mat <- data.matrix(interactions_norm_wide)
        # mat <- mat[, which(colSums(mat) != 0)]
        # mat <- mat[which(rowSums(mat) != 0), ]
        # mat <- mat[order(rownames(mat)), order(colnames(mat))]

        # legend_norm_max <- plyr::round_any(max(mat), 0.01, f = ceiling)
        # create_hm(
        #     mat = mat, col_fun = circlize::colorRamp2(c(0, legend_norm_max), c("white", "red")),
        #     output_file = glue("{args$output_dir}/heatmap_norm__{option}__{current_region}.pdf"),
        #     legend_title = "Fraction detected", save_plot = TRUE,
        #     custom_cell_fun = custom_cell_function_default,
        #     column_title = "Receiver", row_title = "Sender",
        #     column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE
        # )
    }
}
