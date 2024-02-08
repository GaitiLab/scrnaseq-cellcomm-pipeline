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
        description = "Create heatmaps for interactions",
    )
    parser$add_argument("--annot", type = "character", help = "Annotation")
    parser$add_argument("--interactions", type = "character", help = "Interactions")
    parser$add_argument("--agg_level", type = "character", help = "Level of aggregation used: sample or patient")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L1"
    run_name <- "CCI_CellClass_L1"
    args$agg_level <- "sample"
    args$output_dir <- glue("{here::here()}/output/{run_name}/503_heatmaps_interactions")
    args$interactions <- glue("{here::here()}/output/{run_name}/402_aggregation/402_{args$agg_level}_interactions_mvoted_w_filters.rds")
    args$colors <- glue("{here::here()}/000_misc_local/{args$annot}_network_colors.rds")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- glue("{args$output_dir}/{args$agg_level}")
create_dir(output_dir)

# Load additional libraries
pacman::p_load_gh("jokergoo/ComplexHeatmap")

log_info("Load interactions")
interactions <- readRDS(args$interactions)

log_info("Load color dictionary for cell types...")
celltype_colors <- CELLTYPES_COLOR_PALETTE[[args$annot]]


malignant_is_sender <- TRUE

# How to recognize malignant cells
if (args$annot == "CCI_CellClass_L1") {
    malign_id <- "Malignant"
} else {
    malign_id <- "_like"
}

# Titles for heatmap
legend_title <- glue("{str_to_title(args$agg_level)}s")
interaction_axis_title <- "Interaction"

for (malignant_is_sender in c(TRUE, FALSE)) {
    for (option in c("lenient_condition", "lenient_condition_pair")) {
        if (!option %in% colnames(interactions)) {
            next
        }
        # TODO comment later this is only for testing
        # malignant_is_send <- TRUE
        # option <- INTERACTIONS_POST_FILTERING_OPTIONS[1]
        # Set variables
        suffix <- ifelse(malignant_is_sender, "malignant_tme", "tme_malignant")
        malign_function <- ifelse(malignant_is_sender, "source", "target")
        tme_index <- ifelse(malignant_is_sender, 2, 1)
        malign_idx <- ifelse(malignant_is_sender, 1, 2)

        # Pre-filtering based on option
        interactions_filtered <- interactions %>%
            filter(!!sym(option)) %>%
            separate(source_target, c("source", "target"), sep = "__", remove = FALSE) %>%
            mutate(source = factor(source), target = factor(target)) %>%
            rowwise() %>%
            filter(
                # Only interactions between Malignant - TME
                str_detect(!!sym(malign_function), malign_id),
                # No malignant-malignant interactions
                !(str_detect(source, malign_id) && str_detect(target, malign_id))
            )

        # TODO Currently only lenient-voting (based on methods) is implemented
        if (args$agg_level == "sample") {
            interactions_filtered <- interactions_filtered %>%
                mutate(
                    n_detected = str_count(lenient_voting_samples, ",") + 1
                ) %>%
                # Make sure only interactions that have at least 2 samples are included
                filter(n_detected >= 2) %>%
                select(complex_interaction, Region_Grouped, source_target, !!sym(malign_function), n_detected)
            head(interactions_filtered)
            # # A tibble: 6 × 6
            # # Rowwise:
            #   complex_interaction Region_Grouped source          target    source_target              n_detected
            #   <chr>               <fct>          <chr>           <chr>     <chr>                           <dbl>
            # 1 A2M__LRP1           PT             Microglia       Malignant Microglia__Malignant                6
            # 2 A2M__LRP1           TE             Microglia       Malignant Microglia__Malignant                2
            # 3 ADAM10__NRCAM       SC             Oligodendrocyte Malignant Oligodendrocyte__Malignant          3
            # 4 ADAM10__NRCAM       TE             Oligodendrocyte Malignant Oligodendrocyte__Malignant          2
            # 5 AFDN__EPHA7         PT             Malignant       Neuron    Malignant__Neuron                   4
            # 6 AFDN__NRXN3         TE             Malignant       Neuron    Malignant__Neuron                   2
        } else if (args$agg_level == "patient") {
            interactions_filtered <- interactions_filtered %>%
                # Make sure only interactions that have at least 2 samples or patients are included
                filter(lenient_condition_n_patients >= 2) %>%
                select(complex_interaction, Region_Grouped, source_target, !!sym(malign_function), lenient_condition_n_patients) %>%
                rename(n_detected = lenient_condition_n_patients)
        }

        log_info(interactions_filtered %>% pull(complex_interaction) %>% unique() %>% length())

        # Preparing
        df <- interactions_filtered %>%
            mutate(Region_Grouped = factor(Region_Grouped, levels = REGION_GROUPED_LEVELS), source_target = factor(source_target), ) %>%
            arrange(Region_Grouped, !!sym(malign_function), source_target) %>%
            ungroup() %>%
            select(-!!sym(malign_function)) %>%
            pivot_wider(
                names_from = c("Region_Grouped", "source_target", ),
                values_from = n_detected, names_sep = ":", values_fill = 0
            )
        mat <- df %>%
            column_to_rownames("complex_interaction") %>%
            data.matrix()
        annot_regions <- factor(str_split(colnames(mat), ":", simplify = TRUE)[, 1],
            levels = REGION_GROUPED_LEVELS
        )
        # Annotation for cell types
        celltypes <- str_split(colnames(mat), ":", simplify = TRUE)[, 2]
        annot_celltypes <- str_split(celltypes, "__", simplify = TRUE)[, tme_index]
        annot_malign_subtypes <- str_split(celltypes, "__", simplify = TRUE)[, malign_idx]

        # Split into: region - cell type pair
        colnames(mat) <- str_split(colnames(mat), ":", simplify = TRUE)[, 2]
        # split source_target into (1) source - (2) target, keep either 1 or 2
        colnames(mat) <- str_split(colnames(mat), "__", simplify = TRUE)[, tme_index]
        # Use abbreviations
        colnames(mat) <- str_replace_all(colnames(mat), CELLTYPE_ANNOT_ABBREV_DICT)
        # Replace "__" with " - " in interactions (only for visualization purposes)
        rownames(mat) <- str_replace_all(rownames(mat), "__", " - ")

        # Determine max. value of legend
        max_vote <- max(mat)
        col_fun <- circlize::colorRamp2(c(0, 1, max_vote), c("white", "yellow", "blue"))
        cell_size <- 5

        # ---- HEATMAP: INTERACTION X REGION - CELL TYPE PAIR ----
        if (args$annot == "CCI_CellClass_L1") {
            if (malign_function == "source") {
                top_annot <- HeatmapAnnotation(receiver = annot_celltypes, col = list(receiver = celltype_colors))
            } else if (malign_function == "target") {
                top_annot <- HeatmapAnnotation(sender = annot_celltypes, col = list(sender = celltype_colors))
            }
        } else {
            if (malign_function == "source") {
                top_annot <- HeatmapAnnotation(
                    sender = annot_malign_subtypes,
                    receiver = annot_celltypes,
                    col = list(
                        receiver = celltype_colors[names(celltype_colors) %in% annot_celltypes],
                        sender = celltype_colors[names(celltype_colors) %in% annot_malign_subtypes]
                    )
                )
            } else if (malign_function == "target") {
                top_annot <- HeatmapAnnotation(
                    sender = annot_celltypes,
                    receiver = annot_malign_subtypes,
                    col = list(
                        sender = celltype_colors[names(celltype_colors) %in% annot_celltypes],
                        receiver = celltype_colors[names(celltype_colors) %in% annot_malign_subtypes]
                    )
                )
            }
        }
        hm <- Heatmap(mat,
            # Annotation
            column_split = annot_regions,
            top_annotation = top_annot,

            # Clustering
            cluster_column_slices = FALSE,
            cluster_columns = FALSE,
            cluster_rows = TRUE,
            column_gap = unit(2, "mm"),

            # Colors + design
            col = col_fun,
            cell_fun = custom_cell_function_default,

            # Font sizes
            row_names_gp = gpar(fontsize = 12),
            column_names_gp = gpar(fontsize = 12),
            column_title_gp = gpar(fontsize = 12),
            column_names_rot = 45,

            # Size of cells (use square)
            height = nrow(mat) * unit(cell_size, "mm"),
            width = ncol(mat) * unit(cell_size, "mm"),

            # Titles
            name = legend_title,
            row_title = interaction_axis_title,
        )

        hm <- draw(hm)
        hm_size <- get_optimal_output_size(hm)
        output_file <- glue("{output_dir}/503_heatmaps_interactions_{option}_{suffix}.pdf")
        pdf(output_file, width = hm_size$width, height = hm_size$height)
        draw(hm)
        dev.off()

        # ---- FLIPPED HEATMAP: REGION - CELL TYPE PAIR X INTERACTION ----
        if (args$annot == "CCI_CellClass_L1") {
            if (malign_function == "source") {
                row_annot <- rowAnnotation(
                    receiver = annot_celltypes,
                    col = list(receiver = celltype_colors),
                    annotation_legend_param = list(
                        receiver = list(direction = "horizontal", nrow = 1)
                    )
                )
            } else if (malign_function == "target") {
                row_annot <- rowAnnotation(
                    sender = annot_celltypes,
                    col = list(sender = celltype_colors), annotation_legend_param = list(
                        sender = list(direction = "horizontal", nrow = 1)
                    )
                )
            }
        } else {
            if (malign_function == "source") {
                row_annot <- rowAnnotation(
                    sender = annot_malign_subtypes,
                    receiver = annot_celltypes,
                    col = list(
                        receiver = celltype_colors[names(celltype_colors) %in% annot_celltypes],
                        sender = celltype_colors[names(celltype_colors) %in% annot_malign_subtypes]
                    ), annotation_legend_param = list(receiver = list(direction = "horizontal", nrow = 1), sender = list(direction = "horizontal", nrow = 1))
                )
            } else if (malign_function == "target") {
                row_annot <- rowAnnotation(
                    sender = annot_celltypes,
                    receiver = annot_malign_subtypes,
                    col = list(
                        sender = celltype_colors[names(celltype_colors) %in% annot_celltypes],
                        receiver = celltype_colors[names(celltype_colors) %in% annot_malign_subtypes]
                    ), annotation_legend_param = list(
                        sender = list(direction = "horizontal", nrow = 1),
                        receiver = list(direction = "horizontal", nrow = 1)
                    )
                )
            }
        }

        mat <- t(mat)
        mat <- mat[, order(colSums(mat), decreasing = TRUE)]
        hm <- Heatmap(mat,
            # Annotation
            row_split = annot_regions,
            # left_annotation = row_annot,

            # Clustering
            cluster_row_slices = FALSE,
            cluster_columns = TRUE,
            cluster_rows = FALSE,
            row_gap = unit(2, "mm"),

            # Colors + design
            col = col_fun,
            cell_fun = custom_cell_function_default,

            # Font sizes
            row_names_gp = gpar(fontsize = 12),
            column_names_gp = gpar(fontsize = 12),
            column_title_gp = gpar(fontsize = 12),

            # Size of cells (use square)
            height = nrow(mat) * unit(cell_size, "mm"),
            width = ncol(mat) * unit(cell_size, "mm"),

            # Titles
            name = legend_title,
            column_title = interaction_axis_title,
            row_title_rot = 0,
            heatmap_legend_param = list(direction = "horizontal"),
        ) + row_annot

        hm <- draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
        hm_size <- get_optimal_output_size(hm)
        output_file <- glue("{output_dir}/503_heatmaps_interactions_{option}_{suffix}_flipped.pdf")
        pdf(output_file, width = hm_size$width, height = hm_size$height)
        draw(hm)
        dev.off()
    }
}
