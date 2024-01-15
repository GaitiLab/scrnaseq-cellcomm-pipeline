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
        description = "Create Network Visualization",
    )
    parser$add_argument("--input_file", type = "character", help = "Input file")
    parser$add_argument("--interactions_db", type = "character", help = "Interactions database")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L2"
    args$output_dir <- glue("{here::here()}/output/{args$annot}/500_networks")
    args$input_file <- glue("{here::here()}/output/{args$annot}/400_consensus/400_samples_interactions_mvoted_w_filters.rds")
    args$interactions_db <- "001_data_local/interactions_db_v2/ref_db.rds"
    args$meta <- glue("{here::here()}/output/{args$annot}/000_data/gbm_regional_study__metadata.rds")
    args$colors <- glue("{here::here()}/000_misc_local/{args$annot}_network_colors.rds")
    args$min_q <- 0.5
    args$min_cells <- 100
    args$min_celltypes <- 3
    args$has_loops <- FALSE
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
pacman::p_load(randomcoloR, igraph, ggplot2)


log_info("Load interactions...")
interactions <- readRDS(args$input_file)

# Load interactions database
log_info("Determine number of possible interactions...")
interactions_db <- readRDS(args$interactions_db)
n_possible_interactions <- nrow(interactions_db)
log_info(glue("Number of possible interactions: {n_possible_interactions}"))

log_info("Load metadata...")
meta <- readRDS(args$meta)

log_info("Determine number of cells per cell type per region...")
cols_oi <- c("Sample", "Region_Grouped", args$annot)

# Check per sample, the cell types that are included (n_cells >= args$min_cells)
included_celltypes_per_sample <- meta %>%
    select(all_of(cols_oi)) %>%
    group_by_at(vars(all_of(cols_oi))) %>%
    # group_by(all_of(cols_oi)) %>%
    reframe(has_enough_cells = n() >= args$min_cells)

# Number of cell types per sample (with at least args$min_cells)
n_celltypes_per_sample <- included_celltypes_per_sample %>%
    filter(has_enough_cells) %>%
    group_by(Sample, Region_Grouped) %>%
    summarise(n_celltypes = n()) %>%
    filter(n_celltypes >= args$min_celltypes)

included_samples <- n_celltypes_per_sample %>%
    pull(Sample) %>%
    unique()

n_cells_per_type <- meta %>%
    filter(Sample %in% included_samples) %>%
    group_by_at(c("Region_Grouped", "Sample", args$annot)) %>%
    # Number of cells per cell type per sample (per region)
    summarise(n = n()) %>%
    # Deetermine proportion of cells per cell type per sample
    mutate(prop = prop.table(n)) %>%
    ungroup() %>%
    group_by_at(c("Region_Grouped", args$annot)) %>%
    # Determine average proportion of cells per cell type per region
    summarise(avg_prop = mean(prop), avg_pct = 100 * mean(prop))

# TODO, comment after generating colors - DONE (DO NOT TOUCH)
# celltypes <- n_cells_per_type[args$annot] %>%
#     unique() %>% pull()
# celltypes <- celltypes[order(celltypes)]
# colors <- distinctColorPalette(length(celltypes))
# names(colors) <- celltypes
# saveRDS(colors, glue("{here::here()}/000_misc_local/{args$annot}_network_colors.rds"))
colors <- readRDS(args$colors)

# Count the number of interactions between cell type groups
# Types of filters:
# - stringent_region (voting method based stringent + take into account only region)
# - stringent_region_pair (voting stringent + take into account both region and presence of pair in sample of that region)
# - lenient_region (take into account only region)
# - lenient_region_pair (take into account both region and presence of pair in sample of that region)
options <- c("stringent_region", "stringent_region_pair", "lenient_region", "lenient_region_pair")
avail_regions <- interactions %>%
    pull(Region_Grouped) %>%
    unique()
for (option in options) {
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

    for (current_region in avail_regions) {
        # current_region <- "PT"
        log_info("Current selection...")
        df_subset <- interactions_filtered %>%
            filter(Region_Grouped == current_region) %>%
            ungroup() %>%
            select(source, target, n)

        if (!args$has_loops) {
            log_info("Remove loops...")
            df_subset <- df_subset %>% filter(source != target)
        }

        # if (args$only_malignant) {
        #     df_subset <- df_subset %>%
        #         rowwise() %>%
        #         filter(source == "Malignant" || target == "Malignant") %>%
        #         ungroup()
        # }

        df_subset_wide <- df_subset %>% pivot_wider(names_from = target, values_from = n)
        print(df_subset_wide)
        log_info("Determine min and max number of interactions for min-max normalization...")
        edge_list <- df_subset %>%
            # TODO: SQRT transformation - make smaller differences larger
            mutate(n = sqrt(n)) %>%
            # mutate(n = n**2) %>%
            # TODO: Min-max normalization within region
            mutate(n_normalized = (n - min(n)) / (max(n) - min(n))) %>%
            # mutate(n_normalized = n_normalized**2) %>%
            # TODO: Robust normalization within region
            # mutate(n_normalized = (n - quantile(n, 0.25)) / (quantile(n, 0.75) - quantile(n, 0.25))) %>%
            # Transform so that min. value is 0
            # mutate(n_normalized = n_normalized + abs(min(n_normalized))) %>%
            select(source, target, n_normalized) %>%
            ungroup()

        log_info("Create weighted adjacency matrix...")
        adj_mat <- edge_list %>%
            pivot_wider(names_from = target, values_from = n_normalized, values_fill = 0) %>%
            column_to_rownames("source") %>%
            as.matrix()
        adj_mat <- adj_mat[order(rownames(adj_mat)), order(colnames(adj_mat))]

        print(adj_mat)
        log_info("Create igraph network...")
        igraph_network <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "directed", weighted = TRUE, diag = args$has_loops)

        log_info("Add node and edge attributes...")
        noc_subset <- n_cells_per_type %>%
            filter(Region_Grouped == current_region)

        edges <- data.frame(igraph::ends(igraph_network, es = igraph::E(igraph_network), names = TRUE))
        vertices <- data.frame(igraph::V(igraph_network)$name)
        colnames(vertices) <- c("vertex")
        colnames(edges) <- c("sender", "receiver")
        colors_df <- data.frame(colors) %>% rownames_to_column("cell_type")

        edge_colors <- edges %>%
            left_join(colors_df, by = c("sender" = "cell_type")) %>%
            select(colors) %>%
            pull()

        vertex_colors <- vertices %>%
            left_join(colors_df, by = c("vertex" = "cell_type")) %>%
            select(colors) %>%
            pull()

        edge_width <- edges %>%
            left_join(edge_list, by = c("sender" = "source", "receiver" = "target")) %>%
            select(n_normalized) %>%
            pull()

        vertex_size <- vertices %>%
            left_join(noc_subset, by = c("vertex" = args$annot)) %>%
            select(avg_pct) %>%
            pull()

        # Normalize between 0 and 1
        total <- sum(vertex_size)
        vertex_size <- 100 * vertex_size / total

        igraph_network <- igraph_network %>%
            igraph::set_vertex_attr("cell_type", value = igraph::V(igraph_network)$name) %>%
            igraph::set_edge_attr("edge_width", value = edge_width) %>%
            igraph::set_edge_attr("edge_color", value = edge_colors) %>%
            igraph::set_vertex_attr("vertex_color", value = vertex_colors) %>%
            igraph::set_vertex_attr("vertex_size", value = vertex_size)

        # Visualize with igraph
        # Adjust colors to emphasize the most abundant interactions/bigger groups
        vertex_color_adj <- sapply(seq_along(vertex_colors), function(i) {
            adjustcolor(vertex_colors[i], alpha.f = vertex_size[i])
        })
        edge_color_adj <- sapply(seq_along(edge_colors), function(i) {
            if (!is.na(edge_width[i])) {
                if (edge_width[i] > quantile(edge_width, args$min_q)) {
                    adjustcolor(edge_colors[i], alpha.f = 1)
                }
                # if (edge_width[i] > mean(edge_width)) {
                #     adjustcolor(edge_colors[i], alpha.f = 1) }
                else {
                    adjustcolor(edge_colors[i], alpha.f = 0.1)
                }
            } else {
                edge_colors[i]
            }
        })

        # plt_hist <- ggplot() +
        #     geom_histogram(aes(x = edge_width), bins = 10) +
        #     custom_theme() +
        #     geom_vline(xintercept = quantile(edge_width, .25), color = "red") +
        #     geom_vline(xintercept = mean(edge_width), color = "blue") +
        #     geom_vline(xintercept = quantile(edge_width, .50), color = "red") +
        #     geom_vline(xintercept = quantile(edge_width, .75), color = "red")
        # ggsave(glue("{args$output_dir}/hist_{args$annot}_{option}_{current_region}.png"), plt_hist, width = 10, height = 10, units = "in", dpi = 300)
        # auto_crop(glue("{args$output_dir}/hist_{args$annot}_{option}_{current_region}.png"))

        log_info("Plot network...")
        graph_layout <- layout_in_circle(igraph_network)
        has_loops <- ifelse(args$has_loops, "with_loops", "without_loops")
        plot_filename <- glue("{args$output_dir}/network_vis_{args$annot}_{option}_{current_region}_{has_loops}.png")
        png(plot_filename, res = 300, width = 10, height = 10, units = "in")
        plot(igraph_network,
            layout = graph_layout,
            edge.color = edge_color_adj,
            vertex.color = vertex_colors,
            vertex.size = vertex_size,
            edge.curved = 0.1,
            edge.width = edge_width * 10,
            vertex.label.family = "Arial Black",
            vertex.label.color = "black",
            rescale = TRUE,
        )
        dev.off()
        auto_crop(plot_filename)

        # log_info("Create ggnetwork...")
        # nw <- ggnetwork(igraph_network, layout = igraph::layout_in_circle(igraph_network), weights = "weight", niter = 1000, arrow.gap = 0)


        # nw <- nw %>% mutate(name = factor(name), edge_color = factor(edge_color))
        # pt <- nw %>%
        #     ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
        #     geom_edges(aes(linewidth = weight, color = name),
        #         curvature = 0.2, arrow = arrow(length = unit(6, "pt"), type = "closed")
        #     ) +
        #     geom_nodes(aes(size = vertex_size, color = name)) +
        #     custom_theme() +
        #     theme_void() +
        #     geom_nodetext_repel(
        #         color = "black", aes(label = name),
        #         fontface = "bold", size = rel(10)
        #     ) +
        #     theme(legend.position = "none") +
        #     scale_radius(range = c(10, 20)) +
        #     scale_linewidth(range = c(1, 10)) +
        #     scale_colour_manual(values = colors)
        # pt
    }
}
