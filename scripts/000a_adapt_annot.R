if (dir.exists("/opt/.renv")) {
    print("running with Docker/Singularity")
    renv::load("/opt/.renv")
}
print(renv::paths$library())

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
        description = "Adapt cell type annotations",
    )
    parser$add_argument("-i", "--input_file",
        type = "character",
        default = NULL, help = "Path to Seurat object"
    )
    parser$add_argument("-s", "--samples_oi",
        default = NULL,
        help = "Path to sample sheet"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$input_file <- glue("{here::here()}/001_data_local/test_seurat_obj_10x.rds")
    args$output_dir <- glue("{here::here()}/000_misc_local")
    args$samples_oi <- glue("{here::here()}/000_misc_local/samples_oi.txt")
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
pacman::p_load(Seurat, readxl)

log_info("Loading Seurat object...")
seurat_obj <- readRDS(args$input_file)

log_info("Add custom annotation...")
metadata <- seurat_obj@meta.data %>% mutate(CellClass_L4 = case_when(
    startsWith(CellClass_L2, "Malignant") ~ "Malignant",
    TRUE ~ CellClass_L2
))

log_info("Adapt metadata in seurat object...")
seurat_obj@meta.data <- metadata

if (!is.null(args$samples_oi) && file.exists(args$samples_oi)) {
    log_info("Only keeps samples of interest for further processing...")
    samples_oi <- read.table(args$samples_oi, header = TRUE)
    seurat_obj <- subset(seurat_obj, subset = Sample %in% samples_oi$Sample)
} else {
    log_info("No sample sheet provided. All samples will be kept.")
}

log_info("Save adapted object...")
saveRDS(seurat_obj, glue("{output_dir}/seurat_annot_adapted.rds"))
log_info("COMPLETED!")
