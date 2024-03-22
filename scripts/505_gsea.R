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
        description = "Get metadata",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/TESTING")
    args$category <- "H"
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
pacman::p_load(fgsea, msigdbr)

obj <- readRDS("output/CCI_CellClass_L4_w_agg/402_aggregation/402_interactions_combi_agg.rds")
refdb <- readRDS("data/interactions_db/ref_db.rds")
head(refdb)

# REACTOME: 'human' and reactome (subset of C2)
msigdbr_df <- msigdbr(species = "human", category = args$category) %>%
    dplyr::select(gs_name, gene_symbol)

# msigdbr_df <- msigdbr(species = "human", category = "C2") %>%
#     dplyr::select(gs_name, gene_symbol)
# msigdbr_df <- msigdbr_df[startsWith(msigdbr_df$gs_name, "REACTOME_"), ]

# msigdbr_df <- msigdbr(species = "human", category = "C2") %>%
#     dplyr::select(gs_name, gene_symbol)
# msigdbr_df <- msigdbr_df[startsWith(msigdbr_df$gs_name, "KEGG_"), ]

# msigdbr_df <- msigdbr(species = "human", category = "C5") %>%
#     dplyr::select(gs_name, gene_symbol)
# msigdbr_df <- msigdbr_df[startsWith(msigdbr_df$gs_name, "GOBP_"), ]


interactions_to_remove <- obj %>%
    filter(source_target == "Neuron__Progenitor_like") %>%
    pull(complex_interaction) %>%
    unique()

interactions_oi <- obj %>%
    filter(source_target == "Neuron__Neuronal OPC-like") %>%
    pull(complex_interaction) %>%
    unique()

ligands_oi <- obj %>%
    filter(source_target == "Neuron__Neuronal OPC-like", !(complex_interaction %in% interactions_to_remove)) %>%
    left_join(refdb) %>%
    separate(simple_interaction, c("ligand", "receptor")) %>%
    mutate(log10p = -log10(pval)) %>%
    select(ligand, log10p) %>%
    arrange(desc(log10p)) %>%
    distinct(ligand, .keep_all = TRUE)
ligand_oi_names <- ligands_oi %>% pull(ligand)
ligand_oi_ranks <- ligands_oi %>% pull(log10p)
names(ligand_oi_ranks) <- ligand_oi_names

ligands_ref <- obj %>%
    filter(source_target == "Neuron__Progenitor_like") %>%
    left_join(refdb) %>%
    separate(simple_interaction, c("ligand", "receptor")) %>%
    mutate(log10p = -log10(pval)) %>%
    select(ligand, log10p) %>%
    arrange(desc(log10p)) %>%
    distinct(ligand, .keep_all = TRUE)


misgdbr_df_list <- lapply(msigdbr_df %>% pull(gs_name) %>% unique(), function(x) {
    msigdbr_df %>%
        filter(gs_name == x) %>%
        pull(gene_symbol)
})
names(misgdbr_df_list) <- msigdbr_df %>%
    pull(gs_name) %>%
    unique()


fgseaRes <- fgseaMultilevel(
    pathways = misgdbr_df_list, stats = ligand_ranks, minSize = 10,
    maxSize = 1000,
    # nperm = 1e5,
    scoreType = "pos"
)


# look at most significant results
# no significant p-adjusted values
sum(fgseaRes[, padj < 0.1])
sig_results <- fgseaRes[fgseaRes$padj < 0.1, ]
sig_table <- sig_results



topPathwaysUp <- fgseaRes[NES > 0][head(order(padj), n = 20), pathway]
topPathwaysDown <- fgseaRes[NES < 0][head(order(padj), n = 20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(misgdbr_df_list[topPathways], ligand_ranks, fgseaRes,
    gseaParam = 0.5
)
sig_table <- subset(sig_table, NES > 2 | NES < (-2))
ggplot(sig_table, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = NES > 2)) +
    coord_flip() +
    labs(
        x = "Pathway", y = "Normalized Enrichment Score",
        title = "Pathways NES from GSEA"
    )
