# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# Set working directory
set_wd()

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
    args$output_dir <- glue("{here::here()}/output/")
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
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)

options(Seurat.object.assay.version = "v4")
pbmc3k.final <- LoadData("pbmc3k", type = "pbmc3k.final")

sample_ids <- paste0("Sample_", seq(6))

# Manipulate metadata
metadata <- pbmc3k.final@meta.data
metadata$Sample <- sample(sample_ids, size = ncol(pbmc3k.final), replace = TRUE)
metadata <- metadata %>%
    data.frame() %>%
    mutate(
        Condition = case_when(
            Sample %in% sample_ids[1:3] ~ "Group_1",
            Sample %in% sample_ids[4:6] ~ "Group_2"
        ),
        Patient = case_when(
            Sample %in% sample_ids[1:2] ~ "Patient_1",
            Sample %in% sample_ids[3:4] ~ "Patient_2",
            Sample %in% sample_ids[5:6] ~ "Patient_3",
        ),
        seurat_annotations = str_replace_all(seurat_annotations, c(" " = "_", "\\+" = ""))
    )

metadata %>%
    data.frame() %>%
    group_by(Sample, seurat_annotations) %>%
    count() %>%
    pivot_wider(names_from = Sample, values_from = n) %>%
    print()

pbmc3k.final <- AddMetaData(pbmc3k.final, metadata = metadata)


saveRDS(pbmc3k.final, file = "/Users/joankant/Desktop/gaitigroup/Users/Joan/scrnaseq-cellcomm/data/example_data.rds")
