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
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types/514_network_vis_igraph")
    args$input_file <- glue("{here::here()}/output/CellClass_L4_min3_types/400_consensus/400c_post_filtering.rds")
    args$interactions_db <- "001_data_local/interactions_db/interactions_ref.rds"
    args$meta <- glue("{here::here()}/001_data_local/seurat_annot_adapted__metadata.rds")
    args$region_oi <- "TE"
    args$stringency <- 0
    args$colors <- glue("{here::here()}/000_misc_local/colors_CellClassL4.rds")
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
pacman::p_load(randomcoloR, igraph, extrafont)
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
# Determine number of cells per cell type per sample
n_cells_per_type <- meta %>%
    # Number of cells per cell type per sample
    # TODO change if annotation level changes
    count(Region, Sample, CellClass_L4) %>%
    group_by(Region, Sample) %>%
    # Determine proportion of cells per cell type per sample
    mutate(prop = prop.table(n)) %>%
    # Determine proportion of cells per cell type per region (averaging)
    # TODO change if annotation level changes
    group_by(Region, CellClass_L4) %>%
    summarise(avg_n_cells = mean(n), avg_prop_cells = mean(prop), avg_pct_cells = 100 * mean(prop)) %>%
    ungroup()

# Count the number of interactions between cell type groups
interactions_filtered <- interactions %>%
    # TODO comment only for SC region (temporary)
    filter(cond_min_samples_region) %>%
    select(is_stringent, Region, source_target, interaction, source, target) %>%
    group_by(is_stringent, Region, source_target, source, target) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    select(is_stringent, Region, source, target, n) %>%
    ungroup()


# TODO, comment after generating colors
# celltypes <- n_cells_per_type %>%
#     # TODO change if annotation level changes
#     pull(CellClass_L4) %>%
#     unique()
# celltypes <- celltypes[order(celltypes)]
# colors <- distinctColorPalette(length(celltypes))
# names(colors) <- celltypes
# saveRDS(colors, glue("{args$output_dir}/503c_number_of_interactions_network_colors.rds"))

colors <- readRDS(args$colors)

log_info("Current selection...")
df_subset <- interactions_filtered %>%
    filter(is_stringent == args$stringency, Region == args$region_oi) %>%
    ungroup() %>%
    select(-is_stringent, -Region)

df_subset_wide <- df_subset %>% pivot_wider(names_from = target, values_from = n)
print(df_subset_wide)

log_info("Determine min and max number of interactions for min-max normalization...")
max_interactions <- max(df_subset$n)
min_interactions <- min(df_subset$n)
log_info(glue("Min interactions: {min_interactions}"))
log_info(glue("Max interactions: {max_interactions}"))

edge_list <- df_subset %>%
    # Min-max normalization within region
    mutate(n_normalized = (n - min_interactions) / (max_interactions - min_interactions)) %>%
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
igraph_network <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "directed", weighted = TRUE, diag = TRUE)

log_info("Add node and edge attributes...")
noc_subset <- n_cells_per_type %>%
    filter(Region == args$region_oi) # column_to_rownames("CellClass_L4") %>%
# select(avg_pct_cells)

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
    left_join(noc_subset, by = c("vertex" = "CellClass_L4")) %>%
    select(avg_pct_cells) %>%
    pull()

igraph_network <- igraph_network %>%
    igraph::set_vertex_attr("cell_type", value = igraph::V(igraph_network)$name) %>%
    igraph::set_edge_attr("edge_width", value = edge_width) %>%
    igraph::set_edge_attr("edge_color", value = edge_colors) %>%
    igraph::set_vertex_attr("vertex_color", value = vertex_colors) %>%
    igraph::set_vertex_attr("vertex_size", value = vertex_size)



# Visualize with igraph
# plot(igraph_network, edge.width = igraph::E(igraph_network)$weight, main = "circle", vertex.size = igraph::V(igraph_network)$n_cells)

adjustcolor(vertex_colors[1], alpha.f = vertex_size[1] / 100)

# Adjust colors to emphasize the most abundant interactions/bigger groups
vertex_color_adj <- sapply(seq_along(vertex_colors), function(i) {
    adjustcolor(vertex_colors[i], alpha.f = vertex_size[i] / 100)
})
edge_color_adj <- sapply(seq_along(edge_colors), function(i) {
    if (!is.na(edge_width[i])) {
        adjustcolor(edge_colors[i], alpha.f = edge_width[i])
    } else {
        edge_colors[i]
    }
})


log_info("Plot network...")
graph_layout <- layout_in_circle(igraph_network)

extrafont::loadfonts()
pdf(glue("{args$output_dir}/network_vis_{args$region_oi}_stringency_{args$stringency}.pdf"))
plot(igraph_network,
    layout = graph_layout,
    edge.color = edge_color_adj,
    vertex.color = vertex_colors,
    edge.curved = 0.1,
    edge.width = edge_width * 10,
    vertex.size = vertex_size,
    # vertex.label.family = "Arial"
    # vertex.label.dist = 2
)
dev.off()
