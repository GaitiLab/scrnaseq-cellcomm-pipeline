#' Obtain optimal size of heatmap or list of heatmaps
#'
#' @param hm heatmap or list of heatmaps
#' @param m margin
#' @return list of height and width
#' @export
#' @importFrom ComplexHeatmap draw component_height component_width
#' @importFrom grid convertHeight
get_optimal_output_size <- function(hm, m = 4) {
    # First draw the heatmap to obtain the size of the plot
    hm_obj <- draw(hm)

    # Obtain height and width and convert to inches
    ht_height <- sum(component_height(hm_obj)) + unit(m, "mm")
    ht_height <- convertHeight(ht_height, "inch", valueOnly = TRUE)

    ht_width <- sum(component_width(hm_obj)) + unit(m, "mm")
    ht_width <- convertHeight(ht_width, "inch", valueOnly = TRUE)
    return(list(height = ht_height, width = ht_width))
}

#' Cell function for ComplexHeatmap
#'
#' @param j column index
#' @param i row index
#' @param x x coordinate
#' @param y y coordinate
#' @param width width of the cell
#' @param height height of the cell
#' @param fill fill color
#' @return NA
#' @export
#' @importFrom grid grid.rect gpar
custom_layer_function <- function(j, i, x, y, width, height, fill) {
    grid.rect(
        x = x, y = y, width = width, height = height,
        gp = gpar(col = "grey", fill = NA, lwd = 0.2)
    )
}
#' Cell function for ComplexHeatmap
#' @param j column index
#' @param i row index
#' @param x x coordinate
#' @param y y coordinate
#' @param width width of the cell
#' @param height height of the cell
#' @param fill fill color
#' @return NA
#' @export
#' @importFrom grid grid.rect gpar grid.text
custom_cell_function <- function(j, i, x, y, width, height, fill) {
    grid.rect(
        x = x, y = y, width = width, height = height,
        gp = gpar(col = "grey", fill = NA, lwd = 0.2)
    )
    grid.text(mat[i, j], x, y, gp = gpar(fontsize = 5))
}

#' Cell function for ComplexHeatmap
#' @param j column index
#' @param i row index
#' @param x x coordinate
#' @param y y coordinate
#' @param width width of the cell
#' @param height height of the cell
#' @param fill fill color
#' @return NA
#' @export
#' @importFrom grid grid.rect gpar grid.text
custom_cell_function_default <- function(j, i, x, y, width, height, fill) {
    grid.rect(
        x = x, y = y, width = width, height = height,
        gp = gpar(col = "grey", fill = NA, lwd = 0.2)
    )
}

#' Plot heatmap + save
#' @param mat matrix
#' @param col_fun color function
#' @param output_file output file name
#' @param width width of the plot
#' @param height height of the plot
#' @param legend_title legend title
#' @param column_title column title
#' @param column_title_rot rotation of column title
#' @param cluster_columns cluster columns
#' @param show_column_dend show column dendrogram
#' @param show_row_dend show row dendrogram
#' @param save_plot save plot
#' @param col color
#' @param top_annotation column heatmap annotation
#' @return hm
#' @export
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom dplyr %>%
#' @importFrom grid unit gpar grid.rect
#' @importFrom plyr .
create_hm <- function(
    mat,
    col_fun = NULL,
    output_file = "heatmap.pdf",
    legend_title = "score",
    row_title = NULL,
    column_title = NULL, column_title_rot = 90, cluster_columns = TRUE,
    show_column_dend = TRUE, show_row_dend = TRUE, save_plot = FALSE, top_annotation = NULL, show_heatmap_legend = TRUE, heatmap_legend_param = list(), annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE, left_annotation = NULL, custom_cell_fun = NULL, custom_layer_fun = NULL, cluster_rows = TRUE, cell_size = 4, heatmap_legend_list = NULL) {
    hm <- mat %>%
        Heatmap(
            matrix = .,
            col = col_fun,

            # Titles
            name = legend_title,
            column_title = column_title,
            row_title = row_title,

            # Orientation titles
            column_title_rot = column_title_rot,

            # Size of cells (use square)
            height = nrow(.) * unit(cell_size, "mm"),
            width = ncol(.) * unit(cell_size, "mm"),

            # Font sizes
            row_names_gp = gpar(fontsize = 12),
            column_names_gp = gpar(fontsize = 12),
            column_title_gp = gpar(fontsize = 12),

            # Dendrograms
            show_row_dend = show_column_dend,
            show_column_dend = show_row_dend,

            # Clustering
            cluster_columns = cluster_columns,
            cluster_rows = cluster_rows,
            layer_fun = custom_layer_fun,
            cell_fun = custom_cell_fun,
            show_heatmap_legend = show_heatmap_legend,

            # Annotations
            top_annotation = top_annotation,
            left_annotation = left_annotation,
            heatmap_legend_param = heatmap_legend_param
        )
    if (save_plot) {
        hm <- draw(hm, annotation_legend_side = annotation_legend_side, heatmap_legend_side = heatmap_legend_side, merge_legend = merge_legend, heatmap_legend_list = heatmap_legend_list)
        hm_size <- get_optimal_output_size(hm)
        pdf(output_file, width = hm_size$width, height = hm_size$height)
        draw(hm)
        dev.off()
    }
    return(hm)
}

#' Pipeline for creating heatmaps (customized for 401_heatmap)
#' @param interactions_mvoted interactions
#' @param pairs pairs of interest
#' @param setname name of the set
#' @param idx index of the column to use for the title
#' @return list of heatmaps
#' @export
#' @importFrom glue glue
#' @importFrom dplyr %>% filter select group_by summarise pull rowwise mutate select
#' @importFrom tibble column_to_rownames
#' @importFrom stringr str_split str_replace_all
pipeline_heatmap <- function(interactions_mvoted, pairs, obj_literature, setname = "TAM-BDM - Malignant", idx = 2) {
    list_of_hm <- NULL
    all_hm_data <- interactions_mvoted %>%
        filter(source_target %in% pairs)

    # Only keep the interactions that are found for the set of source-targets
    interactions_oi <- all_hm_data %>%
        rowwise() %>%
        mutate(dummy = (across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))) %>%
        select(interaction, dummy) %>%
        group_by(interaction) %>%
        summarise(sm = sum(dummy)) %>%
        filter(sm >= 1) %>%
        pull(interaction)

    # Create a heatmap for each pair, with columns the regions
    for (pair in pairs) {
        hm_mat <- all_hm_data %>%
            filter(source_target == pair, interaction %in% interactions_oi) %>%
            column_to_rownames("interaction") %>%
            select(-source_target) %>%
            as.matrix()
        # Binarise
        hm_mat[hm_mat > 1] <- 1

        # Annotate using literature
        mask <- rownames(hm_mat) %in% (obj_literature %>% filter(source_target == setname) %>% pull(interaction))
        rownames(hm_mat)[mask] <- paste0(rownames(hm_mat)[mask], "*")

        log_info(glue("Number of interactions to plot: {nrow(hm_mat)}"))

        # Add heatmap to list
        list_of_hm <- list_of_hm + create_hm(
            mat = hm_mat,
            legend_title = "Identified",
            column_title = str_split(pair, "__", simplify = TRUE)[idx],
            column_title_rot = 90,
            cluster_columns = FALSE,
            show_column_dend = TRUE, show_row_dend = TRUE,
            save_plot = FALSE, col = colors
        )
    }
    return(list_of_hm)
}

#' Create venndiagram
#' @inheritParams  VennDiagram::venn.diagram
#' @return NA
#' @export
#' @importFrom VennDiagram venn.diagram
#' @importFrom glue glue
#' @importFrom colorjam rainbowJam
create_venndiagram <- function(x, category.names, filename, main = "", ...) {
    venn.diagram(
        x = x,
        category.names = category.names,
        filename = filename,
        output = FALSE,
        main = main,
        disable.logging = FALSE,
        # Output features
        imagetype = "png",
        height = 480,
        width = 480,
        resolution = 600,
        compression = "lzw",

        # Circles
        lwd = 2,
        lty = "blank",
        fill = rainbowJam(length(category.names)),

        # Numbers
        cex = .4,
        fontface = "bold",
        fontfamily = "sans",

        # Set names
        cat.cex = 0.15,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.fontfamily = "sans",

        # Title
        main.fontfamily = "sans",
        main.cex = 0.3,
        main.fontface = "bold"
    )
}

#' Create lineplot of average expression for ligand-receptor pairs
#' @param interactions_df data.frame
#' @return ggplot
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_text scale_alpha_discrete scale_color_manual scale_fill_gradient2 scale_x_discrete scale_y_reverse theme guides
#' @importFrom lemon facet_rep_wrap
create_lineplot <- function(interactions_df) {
    return(ggplot(data = interactions_df, aes(x = type, y = gene_id, group = pair_id)) +
        geom_line(aes(
            x = type, y = gene_id, group = pair_id, color = is_highly_expressed,
            alpha = is_highly_expressed
        ), linewidth = 1) +
        geom_point(aes(x = type, y = gene_id, fill = avg_expr, size = avg_expr * 2), color = "black", pch = 21) +
        # Gene labels
        geom_text(
            data = interactions_df,
            hjust = interactions_df$hjust,
            nudge_x = interactions_df$nudge_x, aes(x = type, y = gene_id, label = genename),
            size = 5
        ) +
        scale_alpha_discrete(range = c(0.1, 1)) +
        scale_color_manual(values = c("yes" = "#A7D489", "no" = "black")) +
        scale_fill_gradient2(
            midpoint = 0, mid = "white", low = "blue", high = "red",
            limits = c(-3, 3), breaks = pretty_breaks()
        ) +
        scale_x_discrete(position = "top") +
        scale_y_reverse() +
        theme(
            axis.line = element_blank(),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), plot.background = element_blank(),
            strip.text = element_blank(),
            legend.title = element_text(angle = 90),
            legend.key.height = unit(1.5, "cm"),
            legend.position = "left"
        ) +
        guides(fill = guide_colorbar(
            title = "Scaled mean expression",
            title.position = "left", title.hjust = 0.5, direction = "vertical",
        ), size = FALSE, colour = FALSE, alpha = FALSE) +
        lemon::facet_rep_wrap(setname ~ source_target, scales = "free_x", repeat.tick.labels = "top"))
}

#' Create lineplot of average expression for ligand-receptor pairs
#' @param avg_expr data.frame
#' @param annot data.frame
#' @param args list of arguments
#' @return list of ggplots
#' @export
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggforce facet_zoom
#' @importFrom dplyr mutate group_by reframe %>%
plot_expr_hist_ccis <- function(avg_expr, annot, args) {
    p <- ggplot(avg_expr, aes(x = expr)) +
        geom_histogram(binwidth = args$binwidth)
    p_data <- ggplot_build(p)$data[[1]]

    matches <- findInterval(annot$expr, p_data$x)
    annot <- annot %>%
        mutate(x_pos = p_data$x[matches]) %>%
        group_by(x_pos) %>%
        reframe(label = paste(gene, collapse = ", "))

    # custom_breaks <- pretty(range(avg_expr$expr),
    #     n = nclass.Sturges(avg_expr$expr),
    #     min.n = 1
    # )

    p_hist <- ggplot(data = avg_expr) +
        geom_histogram(data = avg_expr, aes(x = expr), binwidth = args$binwidth, alpha = 0.5) +
        geom_vline(data = annot, aes(xintercept = x_pos), colour = "red", size = 0.4, linetype = "dashed") +
        geom_text_repel(data = annot, aes(x = x_pos, y = max(p_data$y), label = label), angle = 90, size = 3, hjust = 0.95) +
        custom_theme() +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        facet_zoom(x = expr < quantile(expr, args$max_q))
    p_hist_log_transformed <- p_hist + scale_y_continuous(trans = "log10")
    return(
        list(p_hist, p_hist_log_transformed)
    )
}
