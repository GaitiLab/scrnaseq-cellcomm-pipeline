#' CONSTANT VARIABLES USED IN MULTIPLE SCRIPTS

#' @title List of filtering options
#' @description Types of filtering
#' Count the number of interactions between cell type groups
#' @export
INTERACTIONS_POST_FILTERING_OPTIONS <- c("stringent_condition", "lenient_condition")

# TODO make sure this is up-to-date with metadata conventions in Excel File in OneDrive
#' @title Dictionary of cell type abbreviations
#' @description Abbreviation for the different cell types to be used in visualization
#' @export
CELLTYPE_ANNOT_ABBREV_DICT <- c(
    "Macrophage" = "M",
    "Microglia" = "MG",
    "Oligodendrocyte" = "Oligo",
    "Astrocyte" = "Astro",
    "Neuron" = "N",
    "OPC" = "OPC",
    "Endothelial" = "Endo",
    "Malignant" = "Mal",
    "Progenitor_like" = "PL",
    "Differentiated_like" = "DL",
    "T_cell" = "T cell",
    "Myeloid_Immunosuppressive" = "M-Immuno",
    "Myeloid_Inflammatory" = "M-Inflam",
    "Neuronal OPC-like" = "N-OPC-like"
)

#' @title List of regions in GBM data
#' @description Levels for the Region_grouped variable (ordered)
#' @export
REGION_GROUPED_LEVELS <- c("PT", "TE", "SC")
