#' @title Adapt transparency
#' @param ix index
#' @param groups list of cell types of interest
#' @param alpha_new value for new transparency (0 = transparent, 1 = full color)
adapt_transparency <- function(ix, groups, colors, alpha_new = 0.25) {
    if (names(colors)[ix] %in% groups) {
        return(adjustcolor(colors[ix], alpha.f = 1))
    } else {
        return(adjustcolor(colors[ix], alpha.f = alpha_new))
    }
}

# Cell type labels for CCI_CelLClass L1 and L2
#' @title L1 cell type labels
#' @export
L1 <- c(
    # TME
    "Microglia",
    "Macrophage",
    "Neuron",
    "Pericyte",
    "Oligodendrocyte",
    "OPC",
    "Astrocyte",
    "T_cell",
    # Malignant
    "Malignant"
)
#' @title L2 cell type labels
#' @export
L2 <- c(
    # TME
    "Myeloid_Immunosuppressive",
    "Myeloid_Inflammatory",
    "Neuron",
    "Pericyte",
    "Oligodendrocyte",
    "OPC",
    "Astrocyte",
    "T_cell",
    # Malignant
    "Differentiated_like",
    "Progenitor_like"
)

#' @title L4 cell type labels
#' @export
L4 <- c(
    # Malignant
    "Differentiated_like",
    "Progenitor_like",
    "Neuronal OPC-like",
    # TME
    "Myeloid_Immunosuppressive",
    "Myeloid_Inflammatory",
    "Neuron",
    "Pericyte",
    "Oligodendrocyte",
    "OPC",
    "Astrocyte",
    "T_cell"
)


#' @title L2 TME cell type labels (TME only, no tumor)
#' @export
L2_TME <- c(
    "Myeloid_Immunosuppressive",
    "Myeloid_Inflammatory",
    "Neuron",
    "Pericyte",
    "Oligodendrocyte",
    "OPC",
    "Astrocyte",
    "T_cell"
)

#' @title L1 TME cell type labels (TME only, no tumor)
#' @export
L1_TME <- c(
    "Microglia",
    "Macrophage",
    "Neuron",
    "Pericyte",
    "Oligodendrocyte",
    "OPC",
    "Astrocyte",
    "T_cell"
)
palette <- palette.colors(palette = "Classic Tableau")
# scales::show_col(palette)

#' @title Color palette for L1 TME cell types
#' @export
COLORS_L1_TME <- palette[seq_along(L1_TME)]
names(COLORS_L1_TME) <- L1_TME

#' @title Color palette for L2 TME cell types
#' @export
COLORS_L2_TME <- palette[seq_along(L2_TME)]
names(COLORS_L2_TME) <- L2_TME

# Current groups of interest for highlighting
groups_oi <- c("Neuron")

#' @title Color palette for L1 TME cell types (highlight groups of interest)
#' @export
COLORS_L1_TME_HIGHLIGHTED <- sapply(seq_along(COLORS_L1_TME), adapt_transparency, groups = groups_oi, colors = COLORS_L1_TME)
names(COLORS_L1_TME_HIGHLIGHTED) <- L1_TME

# scales::show_col(COLORS_L1_TME)
# scales::show_col(COLORS_L1_TME_HIGHLIGHTED)

#' @title Color palette for L2 TME cell types (highlight groups of interest)
#' @export
COLORS_L2_TME_HIGHLIGHTED <- sapply(seq_along(COLORS_L2_TME), adapt_transparency, groups = groups_oi, colors = COLORS_L2_TME)
names(COLORS_L2_TME_HIGHLIGHTED) <- L2_TME
# scales::show_col(COLORS_L2_TME)
# scales::show_col(COLORS_L2_TME_HIGHLIGHTED)

# #' @title Color palette for L2 TME cell types (highlight groups of interest)
# #' @export
# COLORS_L4_TME_HIGHLIGHTED <- sapply(seq_along(COLORS_L4_TME), adapt_transparency, groups = groups_oi, colors = COLORS_L4_TME)
# names(COLORS_L4_TME_HIGHLIGHTED) <- L4_TME

#' @title Color palette for L1 cell types
#' @export
COLORS_L1 <- palette[seq_along(L1)]
names(COLORS_L1) <- L1

#' @title Color palette for L2 cell types
#' @export
COLORS_L2 <- palette[seq_along(L2)]
names(COLORS_L2) <- L2


#' @title Color palette for L2 cell types
#' @export
COLORS_L4 <- palette[seq_along(L4)]
names(COLORS_L4) <- L4

# scales::show_col(COLORS_L1)
# scales::show_col(COLORS_L2)
#' @title Dictionary of color palettes
#' @export
CELLTYPES_COLOR_PALETTE <- list("CCI_CellClass_L1" = COLORS_L1, "CCI_CellClass_L2" = COLORS_L2, "CCI_CellClass_L2_2" = COLORS_L4)
