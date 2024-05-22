# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# require(GBMutils)
# Set working directory
set_wd()

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
    parser$add_argument("-n", "--sample_var", type = "character", default = "Sample", help = "Name of sample variable, necessary for splitting")
    parser$add_argument("--downsampling_sheet", type = "character", default = "", help = "Path to RDS file with the cell ids for downsampling")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$input_file <- glue("{here::here()}/output/CCI_CellClass_L1_conf_min50/000_data/split_by_Sample/6419_cortex.rds")
    args$output_dir <- glue("{here::here()}/output/test_dowsampling_implementation")
    args$sample_var <- "Sample"
    args$downsampling_sheet <- glue("{here::here()}/output/downsampling_info.rds")
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
options(Seurat.object.assay.version = "v4")

pacman::p_load(Seurat)

# Initialization
downsampling_sheet <- NULL
run_ids <- NULL

log_info("Loading Seurat object...")
seurat_obj <- readRDS(args$input_file)
print(seurat_obj)
obj_list <- SplitObject(seurat_obj, split.by = args$sample_var)
obj_ids <- names(obj_list)

log_info("Load downsampling sheet...")
if (file.exists(args$downsampling_sheet)) {
    downsampling_sheet <- readRDS(args$downsampling_sheet)
    run_ids <- names(downsampling_sheet)
}
# obj_id <- obj_ids[1]
log_info(glue("Split object by {args$sample_var}..."))
for (obj_id in obj_ids) {
    log_info(glue("Sample: {obj_id}..."))
    obj <- obj_list[[obj_id]]

    log_info(glue("Number of cells in total: {ncol(obj)}..."))
    if (!is.null(run_ids)) {
        mask <- str_detect(run_ids, obj_id)
        run_ids_masked <- run_ids[mask]
        for (run_id in run_ids_masked[1:3]) {
            # run_id <- run_ids_masked[1]
            ix <- str_split(run_id, "__", simplify = TRUE)[3]
            log_info(glue("Downsampling run: {ix}..."))
            obj_subset <- subset(obj, cells = downsampling_sheet[[run_id]])
            log_info(glue("Saving {args$sample_var}: {obj_id}"))
            saveRDS(obj_subset, glue("{args$output_dir}/{run_id}.rds"))
        }
    } else {
        log_info(glue("Saving {args$sample_var}: {obj_id}"))
        saveRDS(
            obj,
            glue("{output_dir}/{obj_id}.rds")
        )
    }
}
log_info("COMPLETED!")
