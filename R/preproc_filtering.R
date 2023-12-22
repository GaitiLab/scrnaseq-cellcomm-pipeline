# ---- User-functions ----
#' Filtering the Seurat object to reduce size of object
#' 1. Only keep genes that are part of the interactions db
#' 2. Only keep cell types with at least N cells (user-defined)
#' @param (Seurat) seurat_obj Seurat object
#' @param (character) annot Annotation to use for filtering
#' @param (integer) min_cells Minimum number of cells per annotation
#' @param (array) genes_oi list of genes of interest from database
#' @return (Seurat) seurat_obj Filtered seurat object
#'
#' @examples seurat_obj <- filtering(seurat_obj, "custom_annot", 300, genes_oi)
#' @export
filtering <- function(seurat_obj, annot, min_cells, genes_oi) {
  cells_annotated <- seurat_obj[[annot]]
  counts_per_label <- table(cells_annotated)
  labels_to_keep <- names(counts_per_label)[counts_per_label >= min_cells]
  cells_to_keep <- rownames(cells_annotated)[cells_annotated[[annot]] %in% labels_to_keep]
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep, features = genes_oi)
  return(seurat_obj)
}
