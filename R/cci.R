#' @title Filtering the Seurat object to reduce size of object
#' @description  1. Only keep genes that are part of the interactions db, 2. Only keep cell types with at least N cells (user-defined)
#' @param (Seurat) seurat_obj Seurat object
#' @param (character) annot Annotation to use for filtering
#' @param (integer) min_cells Minimum number of cells per annotation
#' @param (array) genes_oi list of genes of interest from database
#' @return (Seurat) seurat_obj Filtered seurat object
#' @examples seurat_obj <- filtering(seurat_obj, "custom_annot", 300, genes_oi)
#' @export
filtering <- function(seurat_obj, annot, min_cells) {
    cells_annotated <- seurat_obj[[annot]]
    counts_per_label <- table(cells_annotated)
    labels_to_keep <- names(counts_per_label)[counts_per_label >= min_cells]
    cells_to_keep <- rownames(cells_annotated)[cells_annotated[[annot]] %in% labels_to_keep]
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    return(seurat_obj)
}

#' @title Ranking interactions for consensus using all tools
#' @param cell_type_pair cell type of interest for filtering
#' @param connectome dataframe with connectome results (interactions ordered by descending score)
#' @param sca dataframe with sca results (interactions ordered by descending score)
#' @param cytotalk dataframe with cytotalk results (interactions ordered by descending score)
#' @param natmi dataframe with natmi results (interactions ordered by descending score)
#' @param logfc dataframe with logfc results (interactions ordered by descending score)
#' @param obj_cell2cell dataframe with cell2cell results (interactions ordered by (ascending) p-values)
#' @param obj_cellchat dataframe with cellchat results (interactions ordered by (ascending) p-values)
#' @param obj_cpdb dataframe with cellphonedb results (interactions ordered by (ascending) p-values)
#' @importFrom dplyr filter %>% pull mutate
#' @importFrom RobustRankAggreg rankMatrix aggregateRanks
#' @export
rank_interactions <- function(cell_type_pair, connectome, sca, cytotalk, natmi, logfc, obj_cell2cell, obj_cellchat, obj_cpdb) {
    ordered_interactions <- list(
        connectome = connectome %>%
            filter(source_target == cell_type_pair) %>%
            pull(complex_interaction),
        sca = sca %>%
            filter(source_target == cell_type_pair) %>%
            pull(complex_interaction),
        cytotalk = cytotalk %>%
            filter(source_target == cell_type_pair) %>%
            pull(complex_interaction),
        natmi = natmi %>%
            filter(source_target == cell_type_pair) %>%
            pull(complex_interaction),
        logfc = logfc %>%
            filter(source_target == cell_type_pair) %>%
            pull(complex_interaction),
        cell2cell = obj_cell2cell %>%
            filter(source_target == cell_type_pair) %>%
            pull(complex_interaction),
        cellchat = obj_cellchat %>%
            filter(source_target == cell_type_pair) %>%
            pull(complex_interaction),
        cpdb = obj_cpdb %>%
            filter(source_target == cell_type_pair) %>%
            pull(complex_interaction)
    )

    r <- rankMatrix(ordered_interactions, full = TRUE)
    #                connectome         sca   cytotalk       natmi       logfc cell2cell  cellchat      cpdb
    # LGALS9__PTPRK 0.001605136 0.096308186 0.13643660 0.011235955 0.027287319 0.3103652 0.9557242 0.5319594
    # VSIR__IGSF11  0.003210273 0.104333868 0.08025682 0.009630819 0.014446228        NA 0.9993675 0.9832736
    # TGFB1__ITGB8  0.004815409 0.117174960 0.14125201 0.008025682 0.008025682 0.4568978 0.9658444 0.3416965
    # PDGFB__PDGFRA 0.006420546 0.101123596 0.07062600 0.004815409 0.003210273 0.8428212 0.9917774 0.4447431
    # C3__LRP1      0.008025682 0.006420546 0.13322632 0.165329053 0.006420546 0.7042733 0.9807084 0.1412784
    # SPP1__ITGA9   0.009630819 0.075441413 0.09630819 0.062600321 0.001605136 0.8926512 0.9946237 0.3578256

    ranked_interactions <- aggregateRanks(rmat = r, method = "RRA")
    #                                  Name        Score
    # NRG2__ERBB2:ERBB3   NRG2__ERBB2:ERBB3 2.000351e-07
    # NRG2__ERBB2:ERBB4   NRG2__ERBB2:ERBB4 8.001405e-07
    # ANXA1__FPR2:FPR3     ANXA1__FPR2:FPR3 5.000878e-06
    # THY1__ITGAM:ITGB2   THY1__ITGAM:ITGB2 9.801721e-06
    # THY1__ITGAX:ITGB2   THY1__ITGAX:ITGB2 1.280225e-05
    # ICAM1__ITGAM:ITGB2 ICAM1__ITGAM:ITGB2 1.620285e-05
    return(ranked_interactions %>% mutate(source_target = cell_type_pair))
}

#' @title Combine p-values
#' @description Combine p-values using Fisher's method
#' @param p_values vector of p-values
#' @return combined p-value
#' @export
#' @importFrom metap sumlog
get_p <- function(x) {
    if (length(x) > 1) {
        return(sumlog(x)$p)
    } else {
        return(x)
    }
}
