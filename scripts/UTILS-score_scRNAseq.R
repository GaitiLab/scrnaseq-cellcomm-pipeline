# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
require(GBMutils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Score scRNAseq sample-wise",
    )
    parser$add_argument("--seurat_obj", type = "character", help = "Seurat object")
    parser$add_argument("--sample_id", type = "character", help = "Sample ID")
    parser$add_argument("--sample_varname", type = "character", help = "Column in metadata representing the sample id")
    parser$add_argument("--keep_confident", type = "numeric", help = "Only true for our own data", default = 0)
    parser$add_argument("--annot", type = "character", help = "Annotation column")
    parser$add_argument("--malignant_label", type = "character", help = "Label representing malignant cells, e.g. 'Malignant' in our own GBM data", default = "")
    parser$add_argument("--signatures", type = "character", help = "list of gene signatures (RDS file)")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/score_scRNAseq/our_cohort")
    args$seurat_obj <- "001_data_local/test_seurat_obj_10x.rds"
    args$sample_varname <- "Sample"
    args$keep_confident <- 1
    args$malignant_label <- "Malignant"
    args$annot <- "consensus_cluster_cell_type_merged"
    args$sample_id <- "6509_cortex"
    args$signatures <- "000_misc_local/gene_lists/neftel_signatures.rds"
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
pacman::p_load(Seurat)

log_info("Load scRNAseq Seurat object...")
obj <- readRDS(args$seurat_obj)

log_info(glue("Select sample: {args$sample_id}..."))
obj <- subset(obj, subset = !!sym(args$sample_varname) == args$sample_id)

# Make sure to use SCT for AddModuleScore()
DefaultAssay(obj) <- "SCT"

# Only for our own GBM data
if (!args$keep_confident) {
    log_info("Filtering: keeping confident cells")

    n_confident_cells <- nrow(obj@meta.data %>% data.frame() %>% filter(Confident_Annotation))

    if (n_confident_cells >= min_cells) {
        obj <- subset(obj, subset = Confident_Annotation)
    } else {
        obj <- NULL
    }
}

if (!is.null(obj)) {
    # By default keep all cells
    cells_to_keep <- rownames(obj@meta.data)
    if (args$malignant_label != "") {
        cells_to_keep <- obj@meta.data %>%
            data.frame(!!sym(args$annot) == args$malignant_label) %>%
            rownames()
    }
    if (length(cells_to_keep) > 0) {
        obj <- subset(obj, cells = cells_to_keep)
        log_info("Load signatures...")
        signatures <- readRDS(args$signatures)
        signatures[["invasive_signature"]] <- load_invasive_signature("high")

        log_info("Compute Neftel scores...")
        obj <- AddModuleScore(obj,
            features = signatures,
            name = paste0(names(signatures), "__")
        )

        # Extract scores from metadata + Sample / regional info
        scores <- obj@meta.data %>%
            data.frame() %>%
            select(matches(names(signatures)))

        colnames(scores) <- str_split(colnames(scores), "__", simplify = TRUE)[, 1]

        meta_oi <- obj@meta.data %>%
            data.frame() %>%
            select(Sample, region)

        res <- merge(meta_oi, scores, by = "row.names")

        log_info("Save scores...")
        saveRDS(res, glue("{args$output_dir}/{args$sample_id}_scored.rds"))
    }
}
