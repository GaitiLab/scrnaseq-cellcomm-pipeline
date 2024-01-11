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
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L2"
    args$output_dir <- glue("{here::here()}/output/{args$annot}/501_heatmap")
    args$input_file <- glue("{here::here()}/output/{args$annot}/400_consensus/400_samples_interactions_mvoted_w_filters.rds")
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

log_info("Load interactions...")
interactions <- readRDS(args$input_file)

# Count the number of interactions between cell type groups
# Types of filters:
# - stringent_region (voting method based stringent + take into account only region)
# - stringent_region_pair (voting stringent + take into account both region and presence of pair in sample of that region)
# - lenient_region (take into account only region)
# - lenient_region_pair (take into account both region and presence of pair in sample of that region)
options <- c("stringent_region", "lenient_region", "lenient_region_pair", "stringent_region_pair")
avail_regions <- interactions %>%
    pull(Region_Grouped) %>%
    unique()
for (option in options) {
    for (current_region in avail_regions) {
        # option <- "lenient_region"
        # current_region <- "PT"

        interactions_filtered <- interactions %>%
            ungroup() %>%
            # TODO: filter based on one of the above options
            filter(!!sym(option), !is.na(!!sym(option))) %>%
            select(Region_Grouped, source_target, complex_interaction) %>%
            group_by(Region_Grouped, source_target) %>%
            separate(source_target, c("source", "target"), sep = "__") %>%
            group_by(Region_Grouped, source, target) %>%
            summarise(n = n()) %>%
            ungroup()

        log_info("Current selection...")
        interactions_subset_long <- interactions_filtered %>%
            filter(Region_Grouped == current_region) %>%
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
            output_file = glue("{args$output_dir}/heatmap__{option}__{current_region}.pdf"),
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
