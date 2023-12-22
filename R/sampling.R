#' Sample cells from a Seurat object for a given annotation
#'
#' @param annot An annotation to sample cells from (label, e.g. "Tumor")
#' @param annot_name The name of the annotation column in the Seurat object (e.g. "custom_annot")
#' @param seurat_obj A Seurat object
#' @param n_cells The number of cells to sample for each annotation
#' @return A vector of cell IDs
#' @export
#' @importFrom dplyr %>% pull
sample_cells <- function(annot, annot_name, n_cells, seurat_obj) {
  return(sample(
    colnames(seurat_obj)[seurat_obj[[annot_name]] %>% pull() == annot],
    n_cells
  ))
}

#' Sample cells for an iteration (or run)
#'
#' @param i The iteration number
#' @param annot_name The name of the annotation column in the Seurat object (e.g. "custom_annot")
#' @param cell_types A vector of cell types to sample from
#' @param n_cells The number of cells to sample for each annotation
#' @param seurat_obj A Seurat object
#' @return a square matrix (n_cells x n_cells), columns = cell types, rows = cell IDs
#' @export
sample_iteration <- function(i, annot_name, cell_types, n_cells, seurat_obj) {
  # For each cell type of interest sample the cells.
  return(sapply(cell_types, sample_cells,
    annot_name = annot_name, n_cells = n_cells, seurat_obj = seurat_obj
  ))
}

#' Create a dataframe from a sampling grid
#'
#' @param i The iteration number
#' @param sampling_grid A sampling grid
#' @return A dataframe
#' @export
#' @importFrom dplyr mutate %>%
convert_iteration_matrix_to_df <- function(i, sampling_grid) {
  return(as.data.frame(sampling_grid[[i]]) %>% mutate(run = i))
}
