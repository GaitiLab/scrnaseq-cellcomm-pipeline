# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

pacman::p_load_gh("GaitiLabUtils")

# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Create differential heatmap of number of interactions",
    )
    parser$add_argument("--input_file", type = "character", default = "", help = "Aggregated interactions file named: 402_interactions_combi_agg.rds")
    parser$add_argument("--condition_varname", type = "character", default = "", help = "Column w/ 'condition' with multiple categorical levels to compare, e.g. Mutation, Region")
    parser$add_argument("--group1", type = "character", default = "", help = "Groups to compare: number of interactions in {group1} - {group2}")
    parser$add_argument("--group2", type = "character", default = "", help = "Groups to compare: number of interactions in {group1} - {group2}")
    parser$add_argument("--color_group1", type = "character", default = scales::muted("blue"), help = "Color for group 1 used for color gradient in heatmap; Groups to compare: number of interactions in {group1} - {group2}")
    parser$add_argument("--color_group2", type = "character", default = scales::muted("red"), help = "Color for group 2 used for color gradient in heatmap; Groups to compare: number of interactions in {group1} - {group2}")
    parser$add_argument("--remove_autocrine", type = "numeric", default = 1)
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/TESTING_PIPELINE")
    args$input_file <- "/Users/joankant/Desktop/gaitigroup/Users/Jiaoyi/scrnaseq-cellcomm/output/cci_scvi_merged_annotation_perSample_LP_IMM_Apr9/402_aggregation/402_interactions_combi_agg.rds"
    args$condition_varname <- "Mutation"
    args$group1 <- "BRCA1"
    args$group2 <- "NonCarrier"
    args$remove_autocrine <- 0
    args$color_group2 <- scales::muted("blue")
    args$color_group1 <- scales::muted("red")
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
pacman::p_load(ComplexHeatmap)

log_info("Load data + formatting...")
obj <- readRDS(args$input_file) %>%
    # Additional filtering on p-value (aggregate rank)
    filter(pval < 0.05) %>%
    separate(source_target, c("source", "target"), sep = "__", remove = FALSE) %>%
    # Remove direction by sorting source-target alphabetically
    rowwise() %>%
    mutate(source_target_undirected = paste0(sort(c(source, target)), collapse = "__")) %>%
    # Remove duplicate interactions, when they are found in both directions only keep 1
    distinct(!!sym(args$condition_varname), source_target_undirected, complex_interaction, .keep_all = FALSE) %>%
    separate(source_target_undirected, c("source", "target"), sep = "__") %>%
    arrange(!!sym(args$condition_varname), source, target, .by_group = TRUE) %>%
    # Count number of interactions per source-target per region
    group_by(!!sym(args$condition_varname), target, source) %>%
    summarise(n = n()) %>%
    ungroup()
# head(obj)
# # A tibble: 6 × 4
#   !!sym(args$condition_varname target                    source                     n
#   <fct>          <chr>                     <chr>                  <int>
# 1 PT             Astrocyte                 Astrocyte                101
# 2 PT             Differentiated_like       Differentiated_like        4
# 3 PT             Invasive-high OPC/NPC1    Astrocyte                 48
# 4 PT             Invasive-high OPC/NPC1    Differentiated_like       23
# 5 PT             Invasive-high OPC/NPC1    Invasive-high OPC/NPC1    83
# 6 PT             Myeloid_Immunosuppressive Astrocyte                 92

log_info(glue("Convert to wide format based on levels in '{args$condition_varname}'..."))
obj_wide_condition <- obj %>%
    pivot_wider(names_from = !!sym(args$condition_varname), values_from = n, values_fill = 0) %>%
    select(source, target, !!sym(args$group1), !!sym(args$group2)) %>%
    mutate(diff_n = !!sym(args$group1) - !!sym(args$group2))
# head(obj_wide_condition)
# # A tibble: 6 × 5
#   source                 target                       PT    TC diff_n
#   <chr>                  <chr>                     <int> <int>  <int>
# 1 Astrocyte              Astrocyte                   101     0    101
# 2 Differentiated_like    Differentiated_like           4    77    -73
# 3 Astrocyte              Invasive-high OPC/NPC1       48     0     48
# 4 Differentiated_like    Invasive-high OPC/NPC1       23    67    -44
# 5 Invasive-high OPC/NPC1 Invasive-high OPC/NPC1       83    67     16
# 6 Astrocyte              Myeloid_Immunosuppressive    92     0     92

log_info("Convert to wide format based on cell type labels with the difference in interactions as values...")
obj_diff <- obj_wide_condition %>%
    select(source, target, diff_n) %>%
    arrange(target) %>%
    pivot_wider(names_from = target, values_from = diff_n, values_fill = 0) %>%
    column_to_rownames("source")
# head(obj_diff)
#                           Astrocyte Differentiated_like Invasive-high OPC/NPC1 Myeloid_Immunosuppressive Myeloid_Inflammatory Neuron OPC Oligodendrocyte Progenitor_like T_cell
# Astrocyte                       101                   0                     48                        92                   99    141 125              86              36     15
# Differentiated_like               0                 -73                    -44                      -255                  -33     27   0            -131            -108    -15
# Invasive-high OPC/NPC1            0                   0                     16                       -54                    2    142  42              28             -63     12
# Myeloid_Immunosuppressive         0                   0                      0                       -62                    0     76  62             -62            -194     -1
# Myeloid_Inflammatory              0                   0                      0                         0                   46     80  63              21             -11     26
# Neuron                            0                   0                      0                         0                    0    163 151             133             101     14
log_info("Convert to matrix...")
mat <- data.matrix(obj_diff)

log_info("Make sure order of rows/columns are the same...")
mat <- mat[
    rownames(mat),
    rownames(mat)
]
# print(mat)
#                           Astrocyte Differentiated_like Invasive-high OPC/NPC1 Myeloid_Immunosuppressive Myeloid_Inflammatory Neuron OPC Oligodendrocyte Progenitor_like T_cell
# Astrocyte                       101                   0                     48                        92                   99    141 125              86              36     15
# Differentiated_like               0                 -73                    -44                      -255                  -33     27   0            -131            -108    -15
# Invasive-high OPC/NPC1            0                   0                     16                       -54                    2    142  42              28             -63     12
# Myeloid_Immunosuppressive         0                   0                      0                       -62                    0     76  62             -62            -194     -1
# Myeloid_Inflammatory              0                   0                      0                         0                   46     80  63              21             -11     26
# Neuron                            0                   0                      0                         0                    0    163 151             133             101     14
# OPC                               0                   0                      0                         0                    0      0  61               0              45     19
# Oligodendrocyte                   0                   0                      0                         0                    0      0  84             -31             -83     10
# Progenitor_like                   0                   0                      0                         0                    0      0   0               0             -84     -6
# T_cell                            0                   0                      0                         0                    0      0   0               0               0      3

log_info("Fill whole matrix (mirror)...")
mat <- mat + t(mat)
mat[lower.tri(mat)] <- 0


if (!args$remove_autocrine) {
    log_info("Remove autocrine interactions (diagonal)...")
    # Remove diagonal (autocrine)
    diag(mat) <- 0
    # print(mat)
    #                              Astrocyte Differentiated_like Invasive-high OPC/NPC1 Myeloid_Immunosuppressive Myeloid_Inflammatory Neuron OPC Oligodendrocyte Progenitor_like T_cell
    # Astrocyte                         0                   0                     48                        92                   99    141 125              86              36     15
    # Differentiated_like               0                   0                    -44                      -255                  -33     27   0            -131            -108    -15
    # Invasive-high OPC/NPC1            0                   0                      0                       -54                    2    142  42              28             -63     12
    # Myeloid_Immunosuppressive         0                   0                      0                         0                    0     76  62             -62            -194     -1
    # Myeloid_Inflammatory              0                   0                      0                         0                    0     80  63              21             -11     26
    # Neuron                            0                   0                      0                         0                    0      0 151             133             101     14
    # log_info("Remove fully empty columns/rows...")
} else {
    log_info("Not removing autocrine interactions (diagonal)...")
}
mat <- mat[, which(colSums(mat) != 0)]
mat <- mat[which(rowSums(mat) != 0), ]

# print(mat)
#                           Invasive-high OPC/NPC1 Myeloid_Immunosuppressive Myeloid_Inflammatory Neuron OPC Oligodendrocyte Progenitor_like T_cell
# Astrocyte                                     48                        92                   99    141 125              86              36     15
# Differentiated_like                          -44                      -255                  -33     27   0            -131            -108    -15
# Invasive-high OPC/NPC1                         0                       -54                    2    142  42              28             -63     12
# Myeloid_Immunosuppressive                      0                         0                    0     76  62             -62            -194     -1
# Myeloid_Inflammatory                           0                         0                    0     80  63              21             -11     26
# Neuron                                         0                         0                    0      0 151             133             101     14
# OPC                                            0                         0                    0      0   0              84              45     19
# Oligodendrocyte                                0                         0                    0      0   0               0             -83     10
# Progenitor_like                                0                         0                    0      0   0               0               0     -6


# Determine min/max values for heatmap
legend_max <- plyr::round_any(max(mat), 5, f = ceiling)
legend_min <- plyr::round_any(min(mat), 5, f = floor)

# Setup colors
color_fun <- circlize::colorRamp2(c(legend_min, 0, legend_max), c(args$color_group2, "white", args$color_group1))


chosen_cell_func_w_annot <- get_cell_function(is_upper_tri = TRUE, add_annot = TRUE)

log_info("Plot Heatmap and save...")
create_hm(
    mat = mat,
    col_fun = color_fun,
    output_file = glue("{args$output_dir}/heatmap__diff_undirected_w_annot.pdf"),
    legend_title = glue("{args$group1}-{args$group2}\ninteractions"),
    save_plot = TRUE,
    custom_cell_fun = chosen_cell_func_w_annot,
    column_title = "", row_title = "",
    column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE,
    cell_size = 10
)
chosen_cell_func_w_annot <- get_cell_function(is_upper_tri = TRUE, add_annot = FALSE)

create_hm(
    mat = mat,
    col_fun = color_fun,
    output_file = glue("{args$output_dir}/heatmap__diff_undirected.pdf"),
    legend_title = glue("{args$group1}-{args$group2}\ninteractions"),
    save_plot = TRUE,
    custom_cell_fun = chosen_cell_func_w_annot,
    column_title = "", row_title = "",
    column_title_rot = 0, cluster_rows = FALSE, cluster_columns = FALSE,
    cell_size = 10
)
