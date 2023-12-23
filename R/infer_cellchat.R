
#' Workflow for CellChat
#'
#' @param (Seurat) seurat_obj Seurat object for a sample
#' @param (list) args List of user arguments
#' @param (integer) i Index of sampling grid
#' @return (list) List of CellChat results
#'
#' @export
#' @examples infer_cellchat(seurat_obj, args) or infer_interactions(seurat_obj, args, i = 1)
#' @importFrom glue glue
#' @importFrom reshape2 melt
#' @importFrom pbapply pblapply
#' @importFrom CellChat createCellChat addMeta setIdent subsetData identifyOverExpressedGenes identifyOverExpressedInteractions computeCommunProb filterCommunication aggregateNet computeCommunProbPathway
infer_cellchat <- function(seurat_obj, output_dir, args, i = NULL) {
    log_info("Extract gene expression and convert to matrix...")
    mat <- as.matrix(seurat_obj@assays$RNA@data)
    meta <- seurat_obj@meta.data

    log_info("Create CellChat object...")
    cellchat <- createCellChat(object = mat, meta = meta, group.by = args$ident_col)
    cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = args$ident_col) # set 'labels' as default cell identity

    log_info("Load custom database with interactions...")
    cellchat@DB <- readRDS(args$resource)

    log_info("Preprocessing the expression data...")
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

    future::plan("multisession", workers = args$n_cores)

    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    log_info("Infer cell-cell interactions...")
    cellchat <- computeCommunProb(cellchat, nboot = args$n_perm, population.size = TRUE)
    cellchat <- filterCommunication(cellchat)
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
    saveRDS(res_concat, file = ifelse(is.null(i),
        glue("{output_dir}/cellchat__{get_name(args$gene_expr)}.rds"),
        glue("{output_dir}/cellchat__{get_name(args$gene_expr)}__{i}.rds")
    ))
}
