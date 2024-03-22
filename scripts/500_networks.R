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
    parser$add_argument("--meta", type = "character", help = "Path to metadata")
    parser$add_argument("--min_q", type = "numeric", default = 0.5, help = "Threshold for visibility (transparency)")
    parser$add_argument("--min_cells", type = "numeric", default = 100, help = "Min. cells for a celltype to be included")
    parser$add_argument("--has_loops", type = "numeric", default = 1, help = "Include self-loops")
    parser$add_argument("--annot", type = "character", default = "CCI_CellClass_L1", help = "Annotation variable")
    parser$add_argument("--min_celltypes", type = "numeric", default = 2, help = "Number of min. cell types for a sample to be included")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L2_2"
    run_name <- "CCI_CellClass_L4_w_agg"
    args$output_dir <- glue("{here::here()}/output/{run_name}/500_networks")
    args$input_file <- glue("{here::here()}/output/{run_name}/402_aggregation/402_interactions_combi_agg.rds")
    args$interactions_db <- glue("{here::here()}/data/interactions_db/ref_db.rds")
    args$meta <- glue("{here::here()}/output/{run_name}/000_data/gbm_regional_study__metadata.rds")
    args$min_q <- 0.5
    args$min_cells <- 50
    args$min_celltypes <- 2
    args$has_loops <- TRUE
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- glue("{args$output_dir}")
create_dir(output_dir)

# Load additional libraries
log_info("Load additional libraries...")
pacman::p_load_gh("jokergoo/ComplexHeatmap")
pacman::p_load(igraph, ggplot2)


log_info("Load interactions...")
interactions <- readRDS(args$input_file) %>% filter(pval < 0.05)

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

log_info("Load color dictionary...")

# Setting up colors
colors <- CELLTYPES_COLOR_PALETTE[[args$annot]]

for (option in INTERACTIONS_POST_FILTERING_OPTIONS) {
    if (!option %in% colnames(interactions)) {
        next
    }
    # TODO comment when in use, just for testing
    # option <- INTERACTIONS_POST_FILTERING_OPTIONS[3]
    # Filter on p-value
    interactions_filtered <- interactions %>%
        ungroup() %>%
        # TODO: filter based on one of the above INTERACTIONS_POST_FILTERING_OPTIONS
        filter(!!sym(option), !is.na(!!sym(option))) %>%
        select(Region_Grouped, source_target, complex_interaction) %>%
        group_by(Region_Grouped, source_target) %>%
        separate(source_target, c("source", "target"), sep = "__") %>%
        group_by(Region_Grouped, source, target) %>%
        summarise(n = n()) %>%
        ungroup()
    # head(interactions_filtered)
    #     # A tibble: 6 × 4
    #   Region_Grouped source    target              n
    #   <chr>          <chr>     <chr>           <int>
    # 1 PT             Astrocyte Astrocyte          68
    # 2 PT             Astrocyte Malignant          17
    # 3 PT             Astrocyte Microglia          33
    # 4 PT             Astrocyte Neuron             49
    # 5 PT             Astrocyte OPC                34
    # 6 PT             Astrocyte Oligodendrocyte    22

    for (current_region in REGION_GROUPED_LEVELS) {
        # TODO comment when in use, just for testing
        # current_region <- "PT"
        log_info("Current selection...")
        df_subset <- interactions_filtered %>%
            filter(Region_Grouped == current_region) %>%
            ungroup() %>%
            select(source, target, n)
        if (nrow(df_subset) < 1) {
            next
        }
        # head(df_subset)
        # # A tibble: 6 × 3
        #   source    target              n
        #   <chr>     <chr>           <int>
        # 1 Astrocyte Astrocyte          68
        # 2 Astrocyte Malignant          17
        # 3 Astrocyte Microglia          33
        # 4 Astrocyte Neuron             49
        # 5 Astrocyte OPC                34
        # 6 Astrocyte Oligodendrocyte    22

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
        log_info("Convert to wide format...")
        df_subset_wide <- df_subset %>% pivot_wider(names_from = target, values_from = n)
        # print(df_subset_wide)
        #         # A tibble: 6 × 7
        #   source          Astrocyte Malignant Microglia Neuron   OPC Oligodendrocyte
        #   <chr>               <int>     <int>     <int>  <int> <int>           <int>
        # 1 Astrocyte              68        17        33     49    34              22
        # 2 Malignant              14        89        59     83    26              48
        # 3 Microglia              52        70       151     66    47              51
        # 4 Neuron                 52        69        50    185    58              64
        # 5 OPC                    55        25        18     77    56              39
        # 6 Oligodendrocyte        45        51        62     72    34              50

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
        # r$> head(edge_list)
        # # A tibble: 6 × 3
        #   source    target          n_normalized
        #   <chr>     <chr>                  <dbl>
        # 1 Astrocyte Astrocyte             0.457
        # 2 Astrocyte Malignant             0.0387
        # 3 Astrocyte Microglia             0.203
        # 4 Astrocyte Neuron                0.330
        # 5 Astrocyte OPC                   0.212
        # 6 Astrocyte Oligodendrocyte       0.0962

        log_info("Create weighted adjacency matrix...")
        adj_mat <- edge_list %>%
            pivot_wider(names_from = target, values_from = n_normalized, values_fill = 0) %>%
            column_to_rownames("source") %>%
            as.matrix()
        adj_mat <- adj_mat[order(rownames(adj_mat)), order(colnames(adj_mat))]
        # r$> head(adj_mat)
        #                 Astrocyte  Malignant  Microglia    Neuron Oligodendrocyte       OPC
        # Astrocyte       0.4568600 0.03868717 0.20313826 0.3304670      0.09622478 0.2119000
        # Malignant       0.0000000 0.57732572 0.39955000 0.5445109      0.32318522 0.1376661
        # Microglia       0.3518774 0.46907003 0.86680632 0.4444690      0.34481090 0.3158272
        # Neuron          0.3518774 0.46298712 0.33767480 1.0000000      0.43188878 0.3929198
        # Oligodendrocyte 0.3008725 0.34481090 0.41911043 0.4811069      0.33767480 0.2119000
        # OPC             0.3726786 0.12762337 0.05081063 0.5104871      0.25389331 0.3794856
        log_info("Create igraph network...")
        igraph_network <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "directed", weighted = TRUE, diag = as.logical(args$has_loops))

        log_info("Add node and edge attributes...")
        noc_subset <- n_cells_per_type %>%
            filter(Region_Grouped == current_region)

        edges <- data.frame(igraph::ends(igraph_network, es = igraph::E(igraph_network), names = TRUE))
        vertices <- data.frame(igraph::V(igraph_network)$name)
        colnames(vertices) <- c("vertex")
        colnames(edges) <- c("sender", "receiver")
        colors_df <- data.frame(colors) %>% rownames_to_column("cell_type")
        #          cell_type  colors
        # 1        Astrocyte #D78BC2
        # 2      Endothelial #846FD9
        # 3       Macrophage #D6C9D2
        # 4        Malignant #CDD39C
        # 5        Microglia #D34FD7
        # 6           Neuron #81DCC4
        # 7  Oligodendrocyte #76A9D6
        # 8              OPC #D7D755
        # 9         Pericyte #D87F60
        # 10          T_cell #85E672
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

        log_info("Plot network...")
        graph_layout <- layout_in_circle(igraph_network)
        has_loops <- ifelse(args$has_loops, "with_loops", "without_loops")
        plot_filename <- glue("{output_dir}/network_vis_{args$annot}_{option}_{current_region}_{has_loops}.png")
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
        # auto_crop(plot_filename)

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
