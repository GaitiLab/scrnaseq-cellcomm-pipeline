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
        description = "Create heatmaps",
    )
    parser$add_argument("--input_file", type = "character", help = "Input file")
    parser$add_argument("--interactions_db", type = "character", help = "Interactions database")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types/503c_number_of_interactions")
    args$input_file <- glue("{here::here()}/output/CellClass_L4_min3_types/400_consensus/400c_post_filtering.rds")
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
log_info("Load additional libraries...")
pacman::p_load_gh("jokergoo/ComplexHeatmap")

log_info("Load interactions...")
interactions <- readRDS(args$input_file)

# Load interactions database
log_info("Determine number of possible interactions...")
interactions_db <- readRDS(args$interactions_db)
n_possible_interactions <- nrow(interactions_db)
log_info(glue("Number of possible interactions: {n_possible_interactions}"))

# Count the number of interactions between cell type groups
interactions_filtered <- interactions %>% filter(cond_min_samples_region)
df <- interactions_filtered %>%
    select(is_stringent, Region, source_target, interaction, source, target) %>%
    group_by(is_stringent, Region, source_target, source, target) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    select(is_stringent, Region, source, target, n) %>%
    pivot_wider(names_from = c(target), values_from = n, values_fill = 0)

# For comparison, divide by total possible interactions that can be deteted
df_norm <- df %>% mutate(across(c(interactions_filtered %>% pull(target) %>% unique()), function(x) {
    x / n_possible_interactions
}))

available_regions <- df %>%
    pull(Region) %>%
    unique()
print(available_regions)

stringency_options <- interactions %>%
    pull(is_stringent) %>%
    unique()

# REMARK: only for testing purposes, comment when running
# r <- available_regions[1]
# stringency <- 0

for (stringency in stringency_options) {
    # Determine max value for heatmap legends
    max_fraction <- df_norm %>%
        filter(is_stringent == stringency) %>%
        select(all_of(c(interactions_filtered %>% pull(target) %>% unique()))) %>%
        max()
    legend_norm_max <- plyr::round_any(max_fraction, 0.05, f = ceiling)

    max_n <- df %>%
        filter(is_stringent == stringency) %>%
        select(all_of(c(interactions_filtered %>% pull(target) %>% unique()))) %>%
        max()
    legend_max <- plyr::round_any(max_n, 50, f = ceiling)

    for (r in available_regions) {
        mat <- df %>%
            filter(Region == r, is_stringent == stringency) %>%
            select(-Region, -is_stringent) %>%
            ungroup() %>%
            column_to_rownames("source") %>%
            data.matrix()

        mat <- mat[, which(colSums(mat) != 0)]
        mat <- mat[which(rowSums(mat) != 0), ]
        mat <- mat[order(rownames(mat)), order(colnames(mat))]
        create_hm(
            mat = mat,
            col_fun = circlize::colorRamp2(c(0, legend_max), c("white", "red")), output_file = glue("{args$output_dir}/heatmap__{r}__stringency_{stringency}.pdf"), legend_title = "Interactions",
            save_plot = TRUE,
            custom_cell_fun = custom_cell_function, column_title = "Receiver", row_title = "Sender", column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE
        )

        mat_normalized <- df_norm %>%
            filter(Region == r, is_stringent == stringency) %>%
            select(-Region, -is_stringent) %>%
            ungroup() %>%
            column_to_rownames("source") %>%
            data.matrix()

        # remove rows and columns that are all 0
        mat_normalized <- mat_normalized[, which(colSums(mat_normalized) != 0)]
        mat_normalized <- mat_normalized[which(rowSums(mat_normalized) != 0), ]
        mat_normalized <- mat_normalized[order(rownames(mat_normalized)), order(colnames(mat_normalized))]
        create_hm(mat = mat_normalized, col_fun = circlize::colorRamp2(c(0, legend_norm_max), c("white", "red")), output_file = glue("{args$output_dir}/heatmap__{r}__norm__stringency_{stringency}.pdf"), legend_title = "Fraction detected", save_plot = TRUE, custom_cell_fun = custom_cell_function_default, column_title = "Receiver", row_title = "Sender", column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE)
    }
}
log_info("COMPLETED!")
