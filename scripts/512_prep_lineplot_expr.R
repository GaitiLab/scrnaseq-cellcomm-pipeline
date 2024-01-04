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
        description = "Visualizing average expression of ligands/receptors",
    )
    parser$add_argument("--gene_expr", type = "character", help = "Path to gene expression matrix")
    parser$add_argument("--interactions_ref", type = "character", help = "Path to interactions reference")
    parser$add_argument("--interactions_oi", type = "character", help = "Path to interactions of interest")
    parser$add_argument("--min_pct_expressed", type = "numeric", help = "Minimum percentage of cells expressing a gene", default = 0)
    parser$add_argument("--celltypes_oi", type = "character", help = "Path to celltypes of interest")
    parser$add_argument("--var_annot", type = "character", help = "Cell type variable in object", default = "custom_annot")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$gene_expr <- glue("{here::here()}/final_output/neuron_tumor/100_preprocessing/seurat/6237_2222190_A.rds")
    args$interactions_ref <- glue("{here::here()}/data/interactions_db/interactions_ref.rds")
    args$interactions_oi <- glue("{here::here()}/final_output/neuron_tumor/402_post_filtering/ccis_post_filtering__stringency_0_groupby_1.rds")
    args$output_dir <- glue("{here::here()}/output/neuron_tumor/503_gene_expression")
    args$min_pct_expressed <- 0
    args$celltypes_oi <- glue("{here::here()}/data/celltypes_oi_neuron.txt")

    args$var_annot <- "CellClass_L2"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- paste0(args$output_dir, "/501a_compute_avg_expr")
create_dir(output_dir)

# Additional libraries
pacman::p_load(tidyr, Seurat, dplyr, ggplot2, ggrepel, scales, ggpubr, grid, pbapply)
log_info("Load interactions of interest...")
interactions <- readRDS(args$interactions_oi)

log_info("Loading gene expression...")
gene_expr <- readRDS(args$gene_expr)

region <- gene_expr@meta.data %>%
    pull(region) %>%
    unique()
sample_id <- gene_expr@meta.data %>%
    pull(sample) %>%
    unique()
log_info(glue("region: {region}"))
log_info(glue("sample_id: {sample_id}"))

if (is.null(args$celltypes_oi)) {
    celltypes_oi <- unique(interactions_oi %>% pull(source, target))
} else {
    celltypes_oi <- read.table(args$celltypes_oi, sep = "\t") %>% pull(V1)
}


log_info("Add the main ligand and receptor for each interactions...")
interactions_db <- readRDS(args$interactions_ref) %>% mutate(interaction = str_replace_all(interaction, "__", " - "))
main_lr <- merge(interactions, interactions_db, by = c("interaction"), all.x = TRUE)
main_lr <- main_lr %>% distinct(, .keep_all = TRUE)
log_info(glue("Number of interactions: {nrow(main_lr)}"))

log_info("Compute average expression for all genes in the interactions database...")
all_genes_db <- interactions_db %>%
    select(genename_a, genename_b) %>%
    unlist(use.names = FALSE) %>%
    unique()

expressed_genes_by_type <- pblapply(celltypes_oi, function(celltype) {
    expr_by_type <- subset(gene_expr, subset = (custom_annot == celltype))
    print(expr_by_type)
    if (ncol(expr_by_type) > 0) {
        plot <- Seurat::DotPlot(expr_by_type,
            # Only look at the genes in the interactions database (should already be the case)x
            features = rownames(gene_expr)[rownames(gene_expr) %in% all_genes_db]
        )
        expressed_genes <- plot$data %>%
            dplyr::select(pct.exp, id) %>%
            dplyr::filter(pct.exp > args$min_pct_expressed) %>%
            rownames()
    }
})

names(expressed_genes_by_type) <- celltypes_oi

log_info(glue("Number of interactions after filtering: {nrow(main_lr)}"))

log_info("Compute average expression for all genes in the interactions database...")
out <- lapply(celltypes_oi, function(celltype) {
    gene_expr_avg <- AverageExpression(
        object = subset(gene_expr, custom_annot == celltype), features = expressed_genes_by_type[[celltype]],
        group.by = args$var_annot, slot = "data"
    )$RNA %>%
        as.data.frame() %>%
        rownames_to_column("gene")
    colnames(gene_expr_avg) <- c("gene", celltype)
    return(gene_expr_avg)
})

gene_expr_avg <- out %>%
    reduce(full_join, by = "gene") %>%
    column_to_rownames("gene")


colnames(gene_expr_avg) <- str_replace_all(colnames(gene_expr_avg), "_", " ")
head(gene_expr_avg)

# Standardize (z-scores)
gene_expr_avg_z <- apply(gene_expr_avg, MARGIN = 2, FUN = scale)
rownames(gene_expr_avg_z) <- rownames(gene_expr_avg)

# TODO Hard cut-off.
# msk <- (abs(gene_expr_avg_z) > 3) & !is.na(gene_expr_avg_z)
# gene_expr_avg_z[msk] <- (sign(gene_expr_avg_z[msk]) * 3)

log_debug("Convert wide-format to long-format...")
gene_expr_avg_long <- reshape2::melt(as.matrix(gene_expr_avg_z),
    value.name = "avg_expr", varnames = c("gene", args$var_annot)
)
head(gene_expr_avg_long)
names(expressed_genes_by_type) <- celltypes_oi
gene_expr_avg_long_list <- pblapply(celltypes_oi, function(celltype) {
    return(gene_expr_avg_long %>% filter(custom_annot == celltype, gene %in% expressed_genes_by_type[[celltype]]))
})

# Average expression per gene and cell type
gene_expr_avg_long <- do.call("rbind", gene_expr_avg_long_list)

# log_info("Only keep interactions for which both ligand and receptor are expressed...")
# mask <- sapply(seq_len(nrow(main_lr)), function(i) {
#     ligand <- main_lr[i, "ligand"]
#     receptor <- main_lr[i, "receptor"]
#     source <- main_lr[i, "source"]
#     target <- main_lr[i, "target"]

#     ligand_is_present <- ligand %in% expressed_genes_by_type[[source]]
#     receptor_is_present <- receptor %in% expressed_genes_by_type[[target]]

#     ligand_expr <- nrow(gene_expr_avg_long %>% filter(custom_annot == source, gene == ligand)) > 0

#     receptor_expr <- nrow(gene_expr_avg_long %>% filter(custom_annot == target, gene == receptor)) > 0

#     return(ligand_is_present & receptor_is_present & ligand_expr & receptor_expr)
# })

# main_lr <- main_lr[mask, ]

# log_info("Merge average expression with main ligand-receptor...")
# main_lr$pair_id <- seq_len(nrow(main_lr))
# lr_avg_expr <- merge(main_lr, gene_expr_avg_long,
#     by.x = c("genename_a", "source"), by.y = c("gene", args$var_annot),
#     all.x = TRUE
# )
# head(lr_avg_expr)
# log_debug("Ligand average expression...")
# lr_avg_expr <- lr_avg_expr %>% rename(source_avg_expr_ligand = avg_expr)
# lr_avg_expr <- merge(lr_avg_expr,
#     gene_expr_avg_long,
#     by.x = c("genename_a", "target"),
#     by.y = c("gene", args$var_annot), all.x = TRUE
# )
# lr_avg_expr <- lr_avg_expr %>% rename(target_avg_expr_ligand = avg_expr)

# log_debug("Receptor average expression...")
# lr_avg_expr <- merge(lr_avg_expr, gene_expr_avg_long,
#     by.x = c("genename_b", "source"),
#     by.y = c("gene", args$var_annot), all.x = TRUE
# )
# lr_avg_expr <- lr_avg_expr %>% rename(source_avg_expr_receptor = avg_expr)
# lr_avg_expr <- merge(lr_avg_expr, gene_expr_avg_long,
#     by.x = c("genename_b", "target"), by.y = c("gene", args$var_annot), all.x = TRUE
# )
# lr_avg_expr <- lr_avg_expr %>% rename(target_avg_expr_receptor = avg_expr)
log_info("Save average expression profiles...")
saveRDS(gene_expr_avg_long, glue("{output_dir}/{sample_id}__lr_avg_expr.rds"))

log_info("COMPLETED!")
