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
    args$input_dir <- glue("{here::here()}/final_output/tam-bdm_tumor/201_cci_liana")
    args$output_dir <- glue("{here::here()}/final_output/tam-bdm_tumor")
    # args$samples_oi <- glue("{here::here()}/final_output/tam-bdm_tumor/samples_oi_min200.txt")
    args$celltype_oi <- "Macrophage"
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
pacman::p_load_gh("saezlab/liana")

# TODO Should be replaced with just one metadata object
metadata_parse <- readRDS(glue("{here::here()}/output_tmp/data/metadata__parsebio.rds")) %>% mutate(CellClass_L4 = NA, Dataset = "ParseBio")
metadata_10x <- readRDS(glue("{here::here()}/output_tmp/data/metadata__10x.rds")) %>%
    select(-c(SampleID)) %>%
    mutate(Sample = str_replace_all(Sample, "__", "_"), Dataset = "Multiome 10x")

all_meta <- rbind(metadata_parse, metadata_10x) %>%
    select(Sample, Batch, Patient, Region, Dataset) %>%
    distinct()

# TODO: remove later, as this is temporary for already obtained results
rename_dict <- c("AC-like" = "Malignant_AC", "MES-like" = "Malignant_MES", "NPC & OPC-like" = "Malignant_NPC_&_OPC", "TAM-BDM" = "Macrophage")

log_info("Get filepaths...")
liana_files <- list.files(args$input_dir, full.names = TRUE)

log_info("Load and aggregate all files...")
out <- lapply(liana_files, function(file) {
    obj <- readRDS(file)
    obj_agg <- obj %>%
        liana_aggregate(.) %>%
        mutate(Sample = str_remove(get_name(file), "liana__"))
    return(obj_agg)
})

combined <- do.call(rbind, out)

log_info("Update metadata...")
# Same filtering as in 402_post_filtering.R
combined <- combined %>%
    mutate(source = str_replace_all(source, rename_dict), target = str_replace_all(target, rename_dict)) %>%
    filter(source != target) %>%
    mutate(setname = case_when(
        source == args$celltype_oi ~ glue("{args$celltype_oi} - Malignant"), startsWith(source, "Malignant") & target == args$celltype_oi ~ glue("Malignant - {args$celltype_oi}"), TRUE ~ "Other"
    )) %>%
    filter(setname != "Other")

combined <- combined %>% left_join(all_meta, by = "Sample")

log_info("Save...")
write.csv(combined, glue("{args$output_dir}/505_liana_combined.csv"), row.names = FALSE)

log_info("COMPLETED!")
