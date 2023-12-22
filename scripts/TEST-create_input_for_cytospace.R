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
    args$output_dir <- glue("{here::here()}/output/cytospace_input")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
args <- list()
args$output_dir <- glue("{here::here()}/output/cytospace_input")

create_dir(args$output_dir)

# Load additional libraries
library(Seurat)
library(Matrix)

##########################################################
# README
##########################################################
# Script to generate CytoSPACE input data files
# given the corresponding scRNA Seurat object.
#
# function inputs
#
# 1. scRNA Seurat object
# 2. Path to the output directory to store the results.
#    Default is working directory.
# 3. String to prefix output files with. Default is none.
# 4. Assay to pull counts from. Default is RNA.
##########################################################

generate_cytospace_from_scRNA_seurat_object <- function(scrna_seurat,
                                                        dir_out = ".", fout_prefix = "",
                                                        write_sparse = FALSE, rna_assay = "RNA") {
    scrna_count <- GetAssayData(scrna_seurat, slot = "counts", assay = rna_assay)
    cell_names <- colnames(scrna_count)
    gene_names <- rownames(scrna_count)

    if (dir_out == "") {
        # compatibility with rhe previous version of this script
        message("Note: the output directory will be set to the current working directory.")
    }
    if (!dir.exists(dir_out)) {
        dir.create(dir_out, showWarnings = FALSE)
    }

    if (write_sparse) {
        fout_scrna <- file.path(dir_out, paste0(fout_prefix, "scRNA_data.mtx"))
        fout_genes <- file.path(dir_out, paste0(fout_prefix, "scRNA_data_genes.tsv"))
        fout_cells <- file.path(dir_out, paste0(fout_prefix, "scRNA_data_cells.tsv"))

        Matrix::writeMM(scrna_count, fout_scrna)
        write.table(as.data.frame(gene_names), fout_genes, row.names = F, col.names = F, sep = "\t", quote = F)
        write.table(as.data.frame(cell_names), fout_cells, row.names = F, col.names = F, sep = "\t", quote = F)
    } else {
        scrna_count <- as.data.frame(as.matrix(scrna_count))
        scrna_count <- cbind(rownames(scrna_count), scrna_count)
        colnames(scrna_count)[1] <- "GENES"

        print("Writing output to file")
        fout_scrna <- file.path(dir_out, paste0(fout_prefix, "scRNA_data.txt"))
        write.table(scrna_count, fout_scrna, row.names = F, sep = "\t", quote = F)
    }

    # write cell type labels file
    cell_type_labels <- data.frame(Idents(scrna_seurat))
    rownames(cell_type_labels) <- cell_names
    cell_type_labels <- cbind(rownames(cell_type_labels), cell_type_labels)
    colnames(cell_type_labels) <- c("Cell IDs", "CellType")

    fout_labels <- file.path(dir_out, paste0(fout_prefix, "cell_type_labels.txt"))
    write.table(cell_type_labels, fout_labels, row.names = F, sep = "\t", quote = F)

    print("Done")
}


##########################################################
# README
##########################################################
# Script to generate CytoSPACE input data files
# given the corresponding ST Seurat object.
#
# function inptuts
#
# 1. ST Seurat object
# 2. Path to the output directory to store the results.
#    Default is working directory.
# 3. String to prefix output files with. Default is none.
# 4. Slice name. Default "slice1".
#
##########################################################

generate_cytospace_from_ST_seurat_object <- function(st_seurat,
                                                     dir_out = ".", fout_prefix = "",
                                                     write_sparse = FALSE, slice = "slice1") {
    ST_expressions <- st_seurat@assays$Spatial@counts
    spot_names <- colnames(ST_expressions)
    gene_names <- rownames(ST_expressions)

    if (dir_out == "") {
        # compatibility with the previous version of this script
        message("Note: the output directory will be set to the current working directory.")
    }
    if (!dir.exists(dir_out)) {
        dir.create(dir_out, showWarnings = FALSE)
    }

    if (write_sparse) {
        fout_st <- file.path(dir_out, paste0(fout_prefix, "ST_data.mtx"))
        fout_genes <- file.path(dir_out, paste0(fout_prefix, "ST_data_genes.tsv"))
        fout_spots <- file.path(dir_out, paste0(fout_prefix, "ST_data_cells.tsv"))

        Matrix::writeMM(ST_expressions, fout_st)
        write.table(as.data.frame(gene_names), fout_genes, row.names = F, col.names = F, sep = "\t", quote = F)
        write.table(as.data.frame(spot_names), fout_spots, row.names = F, col.names = F, sep = "\t", quote = F)
    } else {
        ST_expressions <- as.data.frame(as.matrix(ST_expressions))
        ST_expressions <- cbind(rownames(ST_expressions), ST_expressions)
        colnames(ST_expressions)[1] <- "GENES"

        print("Writing output to file")
        fout_st <- file.path(dir_out, paste0(fout_prefix, "ST_data.txt"))
        write.table(ST_expressions, fout_st, row.names = F, sep = "\t", quote = F)
    }

    # write cell type labels file
    coordinates <- st_seurat@images[[slice]]@coordinates[, c("row", "col")]
    coordinates <- cbind(rownames(coordinates), coordinates)
    colnames(coordinates)[1] <- "Spot ID"

    fout_coords <- file.path(dir_out, paste0(fout_prefix, "Coordinates.txt"))
    write.table(coordinates, fout_coords, row.names = F, sep = "\t", quote = F)

    print("Done")
}




# scRNA_Seurat_Object <- readRDS(glue("/cluster/projects/gaitigroup/Users/Yiyan/GBM_10x/multiome_results/Annotated_object/merged.rds"))

# Idents(object = scRNA_Seurat_Object) <- "Malignant.subtyping.L2"
# scRNA_Seurat_Object <- subset(scRNA_Seurat_Object, sample == "6509_enhancing_border")
# generate_cytospace_from_scRNA_seurat_object(scRNA_Seurat_Object, dir_out = args$output_dir, fout_prefix = "6509_enhancing_border", write_sparse = FALSE, rna_assay = "RNA")


obj <- readRDS("/Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-gbm_xenium/data/pimo76_6509_enhancing_border.rds")


obj_df <- data.frame(obj@assays$Xenium@counts) %>% rownames_to_column("V1")

write.table(obj_df, file = "/Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/output/cytospace_input/6509_enhancing_border_STE.txt", sep = "\t")



cells <- fread("/Users/joankant/Desktop/gaitigroup/Data/GBM_Xenium/20230919__202444__Pugh-Shamini-230918/output-XETG00082__0006719__PIMO76-6509-Enhancing_border__20230919__202534/cells.csv.gz")
# generate_cytospace_from_ST_seurat_object(obj, dir_out = glue("{here::here()}/output/cytospace_input"), fout_prefix = "6509_enhancing_border_", write_sparse = FALSE, slice = "slice1")
cells_filtered <- cells[, c("cell_id", "x_centroid", "y_centroid")]
colnames(cells_filtered) <- c("SpotID", "row", "col")

write.table(cells_filtered, file = "/Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/output/cytospace_input/6509_enhancing_border_coordinates.txt")
