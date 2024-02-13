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
        description = "Format annotations for ParseBio",
    )
    parser$add_argument("-i", "--input_file",
        type = "character",
        default = NULL, help = "Path to Seurat object"
    )
    parser$add_argument("--annot", help = "Annotation level to use")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- "5"
    args$input_file <- glue("{here::here()}/output/CCI_CellClass_L1_conf_malign/000_data/split_by_Sample/6234_2895153_A.rds")
    args$output_dir <- glue("{here::here()}/output")
    args$annot <- "CCI_CellClass_L1"
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
pacman::p_load(Seurat)

log_info("Load Seurat object...")
seurat_obj <- readRDS(args$input_file)
DefaultAssay(seurat_obj) <- "RNA"

log_info("Extract RNA assay...")
for (assay_name in names(seurat_obj)) {
    if (assay_name == "RNA") {
        next
    } else {
        log_info(glue("Removing assay: {assay_name}..."))
        try(seurat_obj[[assay_name]] <- NULL, FALSE)
    }
}

# REMARK: Only relevant for GBM project.
# log_info(glue("Number of cells in object, BEFORE: {ncol(seurat_obj)}"))
# seurat_obj <- subset(seurat_obj, subset = Confident_Annotation == TRUE)
# log_info(glue("Number of cells in object, After: {ncol(seurat_obj)}"))

log_info("Saving Seurat object...")
output_name <- str_split(get_name(args$input_file), "__", simplify = TRUE)[1]



saveRDS(seurat_obj, glue("{args$output_dir}/{get_name(output_name)}_reduced_size.rds"))

log_info("COMPLETED!")
