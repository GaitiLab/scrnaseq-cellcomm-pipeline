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
        description = "Combine samples",
    )
    parser$add_argument("--input_dir", required = TRUE, help = "Path to input directory")
    parser$add_argument("--celltype_oi", type = "character", help = "Celltype of interest", default = NULL)
    parser$add_argument("--metadata", type = "character", help = "Path to metadata", default = NULL)
    parser$add_argument("--meta_vars_oi", type = "character", help = "Path to metadata variables of interest", default = NULL)
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$input_dir <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/400_consensus")
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/400_consensus")
    args$celltype_oi <- NULL
    args$metadata <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/000_data/seurat_annot_adapted__metadata.rds")
    args$meta_vars_oi <- ""
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

# log_info("Load interactions - stringent...")
log_info("Load majority voted interactions...")
all_mvoted <- list.files(args$input_dir, full.names = TRUE, pattern = glue("*__interactions_mvoted.rds"))
all_mvoted <- do.call(rbind, lapply(all_mvoted, readRDS))

log_info("Load significant interactions (includes sample information)...")
all_sign_interactions <- list.files(args$input_dir, full.names = TRUE, pattern = glue("*__signif_interactions.rds"))
all_sign_interactions <- do.call(rbind, lapply(all_sign_interactions, readRDS))

log_info("Load metadata...")
if (file.exists(args$meta_vars_oi) && !is.null(args$meta_vars_oi) && args$meta_vars_oi != "") {
    log_info("Metadata variables of interest provided. Using these variables...")
    cols_oi <- read.csv(args$meta_vars_oi, header = FALSE, stringsAsFactors = FALSE) %>%
        pull(V1)
} else {
    log_info("No metadata variables of interest provided. Using default variables (GBM)...")
    cols_oi <- c("Sample", "Patient", "Region", "Batch", "Platform")
}
metadata <- readRDS(args$metadata) %>%
    select(all_of(cols_oi)) %>%
    rownames_to_column("tmp") %>%
    select(-tmp) %>%
    distinct()

log_info("Add metadata...")
all_mvoted <- all_mvoted %>% left_join(metadata, by = "Sample")
all_sign_interactions <- all_sign_interactions %>% left_join(metadata, by = "Sample")

log_info("Save...")
saveRDS(all_mvoted, file = glue("{output_dir}/400_samples_interactions_mvoted.rds"))
saveRDS(all_sign_interactions, file = glue("{output_dir}/400_samples_sign_interactions.rds"))

# if (length(all_files) == 0) {
#     log_error("No files found!")
# }
# combined_samples <- do.call(plyr::rbind.fill, lapply(all_files, function(filepath) {
#     df <- readRDS(filepath)
#     if (str_detect(filepath, "stringency_0")) {
#         df <- df %>%
#             mutate(is_stringent = 0)
#     } else {
#         df <- df %>%
#             mutate(is_stringent = 1)
#     }
#     return(df)
# }))

# all_meta <- readRDS(args$metadata) %>%
#     select(Sample, Batch, Patient, Region) %>%
#     distinct()

# log_info("Add metadata...")
# if (!is.null(args$celltype_oi)) {
#     combined_samples <- combined_samples %>%
#         left_join(all_meta, by = "Sample") %>%
#         select(-source, -target) %>%
#         separate(source_target,
#             into = c("source", "target"), sep = "__", remove = FALSE
#         ) %>%
#         mutate(
#             setname = case_when(
#                 source == args$celltype_oi & startsWith(target, "Malignant") ~ glue("{args$celltype_oi} - Malignant"), startsWith(source, "Malignant") & target == args$celltype_oi ~ glue("Malignant - {args$celltype_oi}"), TRUE ~ "Other"
#             ), interaction = str_replace_all(interaction, "__", " - ")
#         )
# } else {
#     combined_samples <- combined_samples %>%
#         left_join(all_meta, by = "Sample") %>%
#         select(-source, -target) %>%
#         separate(source_target,
#             into = c("source", "target"), sep = "__", remove = FALSE
#         ) %>%
#         mutate(interaction = str_replace_all(interaction, "__", " - "))
# }
# combined_samples <- combined_samples %>%
#     mutate(setname = case_when(str_detect(source, "Malignant") ~ "Malignant-Other", str_detect(target, "Malignant") ~ "Other-Malignant")) %>%
#     mutate(Region = factor(Region, levels = c("PT", "TE", "SC", "NC")))
# log_info("Save...")
# saveRDS(combined_samples, file = glue("{output_dir}/401_samples_combined.rds"))
log_info("COMPLETED!")
