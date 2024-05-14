# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# require(GBMutils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)
if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Inferring CCIs using CellChat"
    )
    parser$add_argument("-p", "--n_perm",
        type = "integer",
        default = 1000,
        help = "Number of permutations for permutation testing"
    )
    parser$add_argument("-db", "--interactions_db",
        type = "character", default = "",
        help = "Path to custom database with interactions (RDS)"
    )
    parser$add_argument("-a", "--annot",
        type = "character", default = "cell_type",
        help = "Column in metadata containing the cell type labels"
    )
    parser$add_argument("-g", "--gene_expr",
        type = "character",
        default = NULL, help = "Seurat object with the gene expression (RDS file)"
    )
    parser$add_argument("-n", "--min_cells",
        type = "integer", default = 5, help = "Minimum number of cells required in each cell group for cell-cell communication"
    )
    parser$add_argument("-nc", "--n_cores",
        default = 1, type = "integer",
        help = "Number of cores to use for parallelization"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/test_pipeline_manual/")
    args$annot <- "seurat_annotations"
    args$n_perm <- 3
    args$interactions_db <- glue("{here::here()}/data/interactions_db/cellchat_db.rds")
    args$gene_expr <- glue("test_pipeline/100_preprocessing/seurat/Sample_2.rds")
    args$n_cores <- 1
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- paste0(args$output_dir)
create_dir(output_dir)

# Load additional libraries
pacman::p_load(reshape2, pbapply)
pacman::p_load_gh("jinworks/CellChat")
options(stringsAsFactors = FALSE)

log_info("Load data...")
seurat_obj <- readRDS(args$gene_expr)
log_info("Extract gene expression and convert to matrix...")
mat <- as.matrix(seurat_obj@assays$RNA@data)
meta <- seurat_obj@meta.data

log_info("Create CellChat object...")
cellchat <- createCellChat(object = mat, meta = meta, group.by = args$annot)
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = args$annot) # set 'labels' as default cell identity

log_info("Load custom database with interactions...")
cellchat@DB <- readRDS(args$interactions_db)

log_info("Preprocessing the expression data...")
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

future::plan("multisession", workers = args$n_cores)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

log_info("Infer cell-cell interactions...")
cellchat <- computeCommunProb(cellchat,
    nboot = args$n_perm,
    population.size = TRUE
)
cellchat <- filterCommunication(cellchat, min.cells = args$min_cells)
# a. signaling pathway level cellchat <- computeCommunProbPathway(cellchat)

# b. aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = glue("{output_dir}/cellchat__{get_name(args$gene_expr)}__raw_obj.rds"))

log_info("Post-processing...")
interactions <- names(cellchat@net$prob[1, 1, ])
res <- pblapply(interactions, function(interaction) {
    # Handle probabilities
    cci <- melt(cellchat@net$prob[, , interaction], )
    colnames(cci) <- c("source", "target", "proba")
    cci["interaction"] <- interaction

    # Handle pvalues
    pval_long <- melt(cellchat@net$pval[, , interaction])
    colnames(pval_long) <- c("source", "target", "pval")
    cci["pval"] <- pval_long$pval
    return(cci)
})
log_info("Concatenate results...")
res_concat <- do.call("rbind", res)

log_info("Save CellChat results...")
saveRDS(
    res_concat,
    glue("{output_dir}/cellchat__{get_name(args$gene_expr)}.rds"),
)

devtools::session_info()
