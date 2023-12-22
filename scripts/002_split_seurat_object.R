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
        description = "Split integrated objects into individual samples (ParseBio)",
    )
    parser$add_argument("-i", "--input_file",
        type = "character",
        default = NULL, help = "Path to input directory"
    )
    parser$add_argument("-n", "--split_varname", type = "character", default = "Sample", help = "Name of sample variable, necessary for splitting")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$input_file <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/001_data/GBM/multiome_10x_reduced_size.rds"
    args$output_dir <- glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/001_data/GBM/samples")
}

# Set up logging
logr <- init_logging(log_level = args$log_level, log_file = "logs/001_prep_split_samples.log")
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- paste0(args$output_dir)
create_dir(output_dir)

# Load additional libraries
options(Seurat.object.assay.version = "v4")

pacman::p_load(Seurat)

log_info("Loading Seurat object...")
seurat_obj <- readRDS(args$input_file)
print(seurat_obj)
obj_list <- SplitObject(seurat_obj, split.by = args$split_varname)

log_info(glue("Split object by {args$split_varname}..."))
for (obj in obj_list) {
    log_info(glue("Saving: {unique(obj@meta.data[args$split_varname])}..."))
    saveRDS(
        obj,
        glue("{output_dir}/{unique(obj@meta.data[args$split_varname])}.rds")
    )
}
log_info("COMPLETED!")
