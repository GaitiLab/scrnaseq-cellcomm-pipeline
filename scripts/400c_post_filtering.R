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
        description = "Post-filtering/formatting",
    )
    parser$add_argument("--input_file",
        type = "character", default = NULL,
        help = "Path to interactions file"
    )
    parser$add_argument("--metadata",
        type = "character", default = NULL,
        help = "Path to metadata file"
    )
    parser$add_argument("--min_cells",
        type = "integer", default = 100,
        help = "Minimum number of cells for a cell type to be retained in a sample"
    )
    parser$add_argument("--min_frac_samples",
        type = "numeric", default = 0.5,
        help = "Minimum fraction of samples for an interaction to be kept"
    )
    parser$add_argument("--annot",
        type = "character", default = "CellClass_L2",
        help = "Annotation to use for filtering"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CCI_CellClass_L2/400_consensus/")
    args$input_file <- glue("{here::here()}/output/CCI_CellClass_L2/400_consensus/400_samples_interactions_mvoted.rds")
    args$min_frac_samples <- 0.5
    args$metadata <- glue("{here::here()}/output/CCI_CellClass_L2/000_data/gbm_regional_study__metadata.rds")
    args$min_cells <- 100
    args$annot <- "CCI_CellClass_L2"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

number_of_interactions_filename <- glue("{args$output_dir}/number_of_interactions.xlsx")
if (file.exists(number_of_interactions_filename)) {
    file.remove(number_of_interactions_filename)
}

# Load additional libraries
pacman::p_load(xlsx, metap)

# TODO might need to change
glob_min_samples <- 1

log_info("Loading input file...")
input_file <- readRDS(args$input_file)

log_info("Loading metadata...")
metadata <- readRDS(args$metadata)
cols_oi <- c("Sample", "Region_Grouped", args$annot)
# Check per sample, the cell types that are included (n_cells >= args$min_cells)
included_celltypes_per_sample <- metadata %>%
    select(all_of(cols_oi)) %>%
    group_by_at(vars(all_of(cols_oi))) %>%
    # group_by(all_of(cols_oi)) %>%
    summarise(n = n() >= args$min_cells)

source_targets_per_sample <- do.call(rbind, lapply(included_celltypes_per_sample %>% pull(Sample) %>% unique(), function(sample_name) {
    data.frame(Sample = sample_name, source_target = apply(expand.grid(
        included_celltypes_per_sample %>% filter(n, Sample == sample_name) %>% pull(as.symbol(args$annot)),
        included_celltypes_per_sample %>% filter(n, Sample == sample_name) %>% pull(as.symbol(args$annot))
    ), 1, paste, collapse = "__"))
}))

# Determine number of available samples per region
n_samples_by_region <- input_file %>%
    ungroup() %>%
    select(Sample, Region_Grouped) %>%
    distinct() %>%
    group_by(Region_Grouped) %>%
    reframe(total_samples_per_region = n()) %>%
    mutate(min_samples_by_region = ceiling(total_samples_per_region * args$min_frac_samples))

# Combine available source-targets + samples per region
n_samples_by_region_pair <- input_file %>%
    ungroup() %>%
    select(Sample, Region_Grouped) %>%
    distinct() %>%
    left_join(source_targets_per_sample, by = "Sample") %>%
    distinct() %>%
    group_by(Region_Grouped, source_target) %>%
    reframe(total_samples_by_region_pair = n()) %>%
    mutate(min_samples_by_region_pair = ceiling(total_samples_by_region_pair * args$min_frac_samples))

threshold_df <- merge(n_samples_by_region, n_samples_by_region_pair, by = "Region_Grouped")

write.xlsx(threshold_df,
    file = number_of_interactions_filename,
    sheetName = "thresholds", append = FALSE
)

input_file_w_thresholds <- input_file %>%
    left_join(threshold_df, by = c("Region_Grouped", "source_target"))

# lenient voting
input_file_lenient <- input_file %>%
    filter(lenient_voting) %>%
    group_by(complex_interaction, Region_Grouped, source_target) %>%
    reframe(
        lenient_voting_n_samples = n(),
        lenient_voting_samples = paste0(Sample, collapse = ", "),
        # pval_combined = get_p(pval)
    )

# stringent voting
input_file_stringent <- input_file %>%
    filter(stringent_voting) %>%
    group_by(complex_interaction, Region_Grouped, source_target) %>%
    reframe(
        stringent_voting_n_samples = n(),
        stringent_voting_samples = paste0(Sample, collapse = ", "),
        # pval_combined = get_p(pval)
    )

# Combine lenient and stringent voting
input_file_voted <- merge(input_file_lenient, input_file_stringent)
cols_to_remove <- c("lenient_voting_n_samples", "stringent_voting_n_samples", "min_samples_by_region", "min_samples_by_region_pair", "total_samples_per_region", "total_samples_by_region_pair")
input_file_w_filters <- merge(input_file_w_thresholds, input_file_voted, all.x = TRUE) %>%
    rowwise() %>%
    mutate(
        # TODO: should we enforce a min. of 2 samples?
        lenient_region = (lenient_voting_n_samples >= min_samples_by_region) && (lenient_voting_n_samples >= glob_min_samples),
        lenient_region_pair = (lenient_voting_n_samples >= min_samples_by_region_pair) && (lenient_voting_n_samples >= glob_min_samples),
        stringent_region = (stringent_voting_n_samples >= min_samples_by_region) && (stringent_voting_n_samples >= glob_min_samples),
        stringent_region_pair = (stringent_voting_n_samples >= min_samples_by_region_pair) && (stringent_voting_n_samples >= glob_min_samples),
        source = str_split(source_target, "__", simplify = TRUE)[, 1],
        target = str_split(source_target, "__", simplify = TRUE)[, 2],
        setname = case_when(
            str_detect(source, "Malignant") ~ "Malignant-Other",
            str_detect(target, "Malignant") ~ "Other-Malignant",
            TRUE ~ "default"
        ),
        Region_Grouped = factor(Region_Grouped, levels = c("PT", "TE", "SC", "NC"))
    ) %>%
    select(-all_of(cols_to_remove))

log_info(glue("Number of interactions: {nrow(input_file)}"))
log_info(glue("Number of interactions: {nrow(input_file_w_filters)}"))

# Save in Excel file

lenient_by_region <- input_file_w_filters %>%
    filter(lenient_region) %>%
    group_by(Region_Grouped, source, target) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = target, values_from = n, values_fill = NA) %>%
    data.frame()
write.xlsx(lenient_by_region,
    file = number_of_interactions_filename,
    sheetName = "lenient_by_region", append = TRUE
)

lenient_by_region_pair <- input_file_w_filters %>%
    filter(lenient_region_pair) %>%
    group_by(Region_Grouped, source, target) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = target, values_from = n, values_fill = NA) %>%
    data.frame()
write.xlsx(lenient_by_region_pair,
    file = number_of_interactions_filename,
    sheetName = "lenient_by_region_pair", append = TRUE
)

stringent_by_region <- input_file_w_filters %>%
    filter(stringent_region) %>%
    group_by(Region_Grouped, source, target) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = target, values_from = n, values_fill = NA) %>%
    data.frame()
write.xlsx(stringent_by_region,
    file = number_of_interactions_filename,
    sheetName = "stringent_by_region", append = TRUE
)

stringent_by_region_pair <- input_file_w_filters %>%
    filter(stringent_region_pair) %>%
    group_by(Region_Grouped, source, target) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = target, values_from = n, values_fill = NA) %>%
    data.frame()
write.xlsx(stringent_by_region_pair,
    file = number_of_interactions_filename,
    sheetName = "stringent_by_region_pair", append = TRUE
)

log_info("Save output...")
saveRDS(input_file_w_filters,
    file = glue("{args$output_dir}/400_samples_interactions_mvoted_w_filters.rds")
)
log_info("COMPLETED!")
