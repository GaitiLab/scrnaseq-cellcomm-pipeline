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
    parser$add_argument("--min_celltypes",
        type = "integer", default = 3,
        help = "Minimum number of cell types for a sample to be retained"
    )
    parser$add_argument("--sample_varname", type = "character", help = "Name of sample variable", default = "Sample")
    parser$add_argument("--condition_varname", type = "character", help = "Name of condition variable", default = "Region_Grouped")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L2"
    args$output_dir <- glue("{here::here()}/output/{args$annot}/400_consensus/")
    args$input_file <- glue("{here::here()}/output/{args$annot}/400_consensus/400_samples_interactions_mvoted.rds")
    args$min_frac_samples <- 0.5
    args$metadata <- glue("{here::here()}/output/{args$annot}/000_data/gbm_regional_study__metadata.rds")
    args$min_cells <- 100
    args$min_celltypes <- 3
    args$sample_varname <- "Sample"
    args$condition_varname <- "Region_Grouped"
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

# TODO might need to change this to too, to not get a skewed view of the data
glob_min_samples <- 1

log_info("Loading input file...")
input_file <- readRDS(args$input_file)
# r$> head(input_file)
# # A tibble: 6 × 11
#   Sample         source_target        complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb lenient_voting stringent_voting Region_Grouped
#   <chr>          <chr>                <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl> <lgl>          <lgl>            <chr>
# 1 6234_2895153_A Microglia__Microglia A2M__LRP1                   4        1           1            1       1 TRUE           TRUE             TE
# 2 6234_2895153_A Microglia__Microglia ACTR2__ADRB2                1        0           0            1       0 FALSE          FALSE            TE
# 3 6234_2895153_A Microglia__Microglia ADAM10__AXL                 4        1           1            1       1 TRUE           TRUE             TE
# 4 6234_2895153_A Microglia__Microglia ADAM10__CADM1               1        0           0            1       0 FALSE          FALSE            TE
# 5 6234_2895153_A Microglia__Microglia ADAM10__CD44                2        0           0            1       1 FALSE          FALSE            TE
# 6 6234_2895153_A Microglia__Microglia ADAM10__GPNMB               4        1           1            1       1 TRUE           TRUE             TE


# COMMENT: Temporary fix, remove later for fully new pipeline runs, correction of interactions database from 7K to 5.4K interactions.
ref_db <- readRDS("001_data_local/interactions_db_v2/ref_db.rds")
input_file <- input_file %>% filter(complex_interaction %in% (ref_db %>% pull(complex_interaction)))

log_info("Loading metadata...")
metadata <- readRDS(args$metadata)
cols_oi <- c(args$sample_varname, args$condition_varname, args$annot)
# Check per sample, the cell types that are included (n_cells >= args$min_cells)
included_celltypes_per_sample <- metadata %>%
    select(all_of(cols_oi)) %>%
    group_by_at(vars(all_of(cols_oi))) %>%
    summarise(has_enough_cells = n() >= args$min_cells) %>%
    ungroup()
# r$> head(included_celltypes_per_sample)
# # A tibble: 6 × 4
#   Sample         Region_Grouped CCI_CellClass_L2    has_enough_cells
#   <chr>          <chr>          <chr>               <lgl>
# 1 6234_2895153_A TE             Astrocyte           FALSE
# 2 6234_2895153_A TE             Differentiated_like FALSE
# 3 6234_2895153_A TE             Endothelial         FALSE
# 4 6234_2895153_A TE             Macrophage          FALSE
# 5 6234_2895153_A TE             Microglia           TRUE
# 6 6234_2895153_A TE             Neuron              TRUE

# Number of cell types per sample (with at least args$min_cells)
cols_oi <- c(args$sample_varname, args$condition_varname)
n_celltypes_per_sample <- included_celltypes_per_sample %>%
    filter(has_enough_cells) %>%
    group_by_at(vars(all_of(cols_oi))) %>%
    summarise(n_celltypes = n()) %>%
    filter(n_celltypes >= args$min_celltypes) %>%
    ungroup()
# r$> head(n_celltypes_per_sample)
# # A tibble: 6 × 3
#   Sample         Region_Grouped n_celltypes
#   <chr>          <chr>                <int>
# 1 6234_2895153_A TE                       3
# 2 6234_2895153_B PT                       4
# 3 6237_2222190_A PT                       7
# 4 6237_2222190_C TE                       5
# 5 6237_2222190_D SC                       6
# 6 6240_4964941_B TE                       3

# Final list of samples to include, should be the same samples as listed in the preprocessing folder.
included_samples <- n_celltypes_per_sample %>%
    pull(Sample) %>%
    unique()
log_info(glue("Number of samples: {length(included_samples)}"))
# r$> head(included_samples)
# [1] "6234_2895153_A" "6234_2895153_B" "6237_2222190_A" "6237_2222190_C" "6237_2222190_D" "6240_4964941_B"

# CHECK: the number of samples in pre-processing should be the same as the number of samples with at least args$min_celltypes cell types
# n_preprocessing <- get_name(list.files("output/CCI_CellClass_L2/100_preprocessing/seurat"))
# log_info(glue("Number of samples in preprocessing: {n_preprocessing %>% length()}"))
# log_info(glue("included_samples == n_preprocessing: {length(intersect(included_samples, n_preprocessing))}"))
# included_celltypes_per_sample <- included_celltypes_per_sample %>% filter(Sample %in% included_samples)
# log_info(glue("Number of samples with at least {args$min_celltypes} cell types: {length(included_samples)}"))

# Potential source-targets based on presence of cell types in individual samples
# Based on the present cell types in each sample, we can determine the potential source-targets
source_targets_per_sample <- do.call(rbind, lapply(included_celltypes_per_sample %>% pull(Sample) %>% unique(), function(sample_name) {
    data.frame(Sample = sample_name, source_target = apply(expand.grid(
        included_celltypes_per_sample %>% filter(has_enough_cells, Sample == sample_name) %>% pull(as.symbol(args$annot)),
        included_celltypes_per_sample %>% filter(has_enough_cells, Sample == sample_name) %>% pull(as.symbol(args$annot))
    ), 1, paste, collapse = "__"))
}))

# Determine number of available samples per region
cols_oi <- c(args$sample_varname, args$condition_varname)
n_samples_by_condition <- input_file %>%
    ungroup() %>%
    select(all_of(cols_oi)) %>%
    distinct() %>%
    group_by_at(args$condition_varname) %>%
    reframe(total_samples_per_condition = n()) %>%
    mutate(min_samples_by_condition = ceiling(total_samples_per_condition * args$min_frac_samples)) %>%
    ungroup()
# r$> head(n_samples_by_condition)
# # A tibble: 3 × 3
#   Region_Grouped total_samples_per_condition min_samples_by_condition
#   <chr>                             <int>                 <dbl>
# 1 PT                                    8                     4
# 2 SC                                   10                     5
# 3 TE                                    8                     4

# Combine available source-targets + samples per region
groupby_cols <- c(args$condition_varname, "source_target")
n_samples_by_condition_pair <- input_file %>%
    ungroup() %>%
    select(all_of(cols_oi)) %>%
    distinct() %>%
    left_join(source_targets_per_sample, by = args$sample_varname) %>%
    distinct() %>%
    group_by_at(vars(all_of(groupby_cols))) %>%
    reframe(total_samples_by_condition_pair = n()) %>%
    mutate(min_samples_by_condition_pair = ceiling(total_samples_by_condition_pair * args$min_frac_samples)) %>%
    ungroup()
# r$> head(n_samples_by_condition_pair)
# # A tibble: 6 × 4
#   Region_Grouped source_target                  total_samples_by_condition_pair min_samples_by_condition_pair
#   <chr>          <chr>                                                 <int>                      <dbl>
# 1 PT             Astrocyte__Astrocyte                                      4                          2
# 2 PT             Astrocyte__Differentiated_like                            1                          1
# 3 PT             Astrocyte__Endothelial                                    1                          1
# 4 PT             Astrocyte__Macrophage                                     1                          1
# 5 PT             Astrocyte__Microglia                                      4                          2
# 6 PT             Astrocyte__Neuron                                         4                          2

# Combine the two tables above,
# total_samples_per_condition and min_samples_by_condition are the same for each source-target
# total_samples_by_condition_pair and min_samples_by_condition_pair are unique for each row (dependent on region)
threshold_df <- merge(n_samples_by_condition, n_samples_by_condition_pair, by = args$condition_varname)
# r$> head(threshold_df)
#   Region_Grouped total_samples_per_condition min_samples_by_condition                  source_target total_samples_by_condition_pair min_samples_by_condition_pair
# 1             PT                        8                     4           Astrocyte__Astrocyte                            4                          2
# 2             PT                        8                     4 Astrocyte__Differentiated_like                            1                          1
# 3             PT                        8                     4         Astrocyte__Endothelial                            1                          1
# 4             PT                        8                     4          Astrocyte__Macrophage                            1                          1
# 5             PT                        8                     4           Astrocyte__Microglia                            4                          2
# 6             PT                        8                     4              Astrocyte__Neuron                            4                          2

log_info("Add to Excel Sheet for future reference...")
write.xlsx(threshold_df,
    file = number_of_interactions_filename,
    sheetName = "thresholds", append = FALSE
)

# ---- (1) Lenient voting ---- #
# See how many samples have identified the same interaction (x region combination)
cols_oi <- c("complex_interaction", args$condition_varname, "source_target")
input_file_lenient <- input_file %>%
    # Only keep the interactions according to the lenient voting
    filter(lenient_voting) %>%
    group_by_at(vars(all_of(cols_oi))) %>%
    reframe(
        lenient_voting_n_samples = n(),
        lenient_voting_samples = paste0(!!sym(args$sample_varname), collapse = ", "),
        # pval_combined = get_p(pval)
    )
# r$> head(input_file_lenient)
# # A tibble: 6 × 5
#   complex_interaction Region_Grouped source_target                    lenient_voting_n_samples lenient_voting_samples
#   <chr>               <chr>          <chr>                                               <int> <chr>
# 1 A2M__LRP1           PT             Endothelial__Astrocyte                                  1 6419_enhancing_border
# 2 A2M__LRP1           PT             Endothelial__Differentiated_like                        1 6419_enhancing_border
# 3 A2M__LRP1           PT             Endothelial__Microglia                                  1 6419_enhancing_border
# 4 A2M__LRP1           PT             Endothelial__Pericyte                                   1 6419_enhancing_border
# 5 A2M__LRP1           PT             Endothelial__Progenitor_like                            1 6419_enhancing_border
# 6 A2M__LRP1           PT             Microglia__Astrocyte                                    3 6237_2222190_A, 6419_cortex, 6509_cortex

# ---- (2) Stringent voting ---- #
input_file_stringent <- input_file %>%
    filter(stringent_voting) %>%
    group_by_at(vars(all_of(cols_oi))) %>%
    reframe(
        stringent_voting_n_samples = n(),
        stringent_voting_samples = paste0(!!sym(args$sample_varname), collapse = ", "),
        # pval_combined = get_p(pval)
    )
# r$> head(input_file_stringent)
# # A tibble: 6 × 5
#   complex_interaction Region_Grouped source_target                    stringent_voting_n_samples stringent_voting_samples
#   <chr>               <chr>          <chr>                                                 <int> <chr>
# 1 A2M__LRP1           PT             Endothelial__Astrocyte                                    1 6419_enhancing_border
# 2 A2M__LRP1           PT             Endothelial__Differentiated_like                          1 6419_enhancing_border
# 3 A2M__LRP1           PT             Endothelial__Microglia                                    1 6419_enhancing_border
# 4 A2M__LRP1           PT             Endothelial__Pericyte                                     1 6419_enhancing_border
# 5 A2M__LRP1           PT             Endothelial__Progenitor_like                              1 6419_enhancing_border
# 6 A2M__LRP1           PT             Microglia__Astrocyte                                      3 6237_2222190_A, 6419_cortex, 6509_cortex

# Combine lenient and stringent voting
# stringent_voting should be a subset of the lenient voting, i.e. number of expected rows is equal to the lenient voting
input_file_voted <- merge(input_file_lenient, input_file_stringent, all = TRUE, by = c("complex_interaction", args$condition_varname, "source_target"))
# head(input_file_voted)
# r$> head(input_file_voted)
#   complex_interaction Region_Grouped                    source_target lenient_voting_n_samples                   lenient_voting_samples stringent_voting_n_samples                 stringent_voting_samples
# 1           A2M__LRP1             PT           Endothelial__Astrocyte                        1                    6419_enhancing_border                          1                    6419_enhancing_border
# 2           A2M__LRP1             PT Endothelial__Differentiated_like                        1                    6419_enhancing_border                          1                    6419_enhancing_border
# 3           A2M__LRP1             PT           Endothelial__Microglia                        1                    6419_enhancing_border                          1                    6419_enhancing_border
# 4           A2M__LRP1             PT            Endothelial__Pericyte                        1                    6419_enhancing_border                          1                    6419_enhancing_border
# 5           A2M__LRP1             PT     Endothelial__Progenitor_like                        1                    6419_enhancing_border                          1                    6419_enhancing_border
# 6           A2M__LRP1             PT             Microglia__Astrocyte                        3 6237_2222190_A, 6419_cortex, 6509_cortex                          3 6237_2222190_A, 6419_cortex, 6509_cortex
# cols_to_remove <- c("lenient_voting_n_samples", "stringent_voting_n_samples", "min_samples_by_condition", "min_samples_by_condition_pair", "total_samples_per_condition", "total_samples_by_condition_pair")

input_file_w_filters <- input_file_voted %>%
    left_join(threshold_df) %>%
    rowwise() %>%
    mutate(
        # TODO: should we enforce a min. of 2 samples?
        lenient_condition = (lenient_voting_n_samples >= min_samples_by_condition) && (lenient_voting_n_samples >= glob_min_samples),
        lenient_condition_pair = (lenient_voting_n_samples >= min_samples_by_condition_pair) && (lenient_voting_n_samples >= glob_min_samples),
        stringent_condition = (stringent_voting_n_samples >= min_samples_by_condition) && (stringent_voting_n_samples >= glob_min_samples),
        stringent_condition_pair = (stringent_voting_n_samples >= min_samples_by_condition_pair) && (stringent_voting_n_samples >= glob_min_samples),
        source = str_split(source_target, "__", simplify = TRUE)[, 1],
        target = str_split(source_target, "__", simplify = TRUE)[, 2],
        # TODO This is only relevant for GBM
        setname = case_when(
            str_detect(source, "Malignant") ~ "Malignant-Other",
            str_detect(target, "Malignant") ~ "Other-Malignant",
            str_detect(source, "_like") ~ "Malignant-Other",
            str_detect(target, "_like") ~ "Other-Malignant",
            TRUE ~ "default"
        ),
    ) %>%
    ungroup()
input_file_w_filters[[args$condition_varname]] <- factor(input_file_w_filters[[args$condition_varname]], levels = c("PT", "TE", "SC", "NC"))

# Check if condition_varname, i.e. region_varname they are ordered
# r$> head(input_file_w_filters[[args$condition_varname]])
# [1] PT PT PT PT PT PT
# Levels: PT TE SC NC

log_info(glue("Number of interactions: {nrow(input_file)}"))
log_info(glue("Number of interactions: {nrow(input_file_w_filters)}"))

# Save in Excel file
groupby_cols <- c(args$condition_varname, "source", "target")
lenient_by_condition <- input_file_w_filters %>%
    filter(lenient_condition) %>%
    group_by_at(vars(all_of(groupby_cols))) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = target, values_from = n, values_fill = NA) %>%
    data.frame()
write.xlsx(lenient_by_condition,
    file = number_of_interactions_filename,
    sheetName = "lenient_by_condition", append = TRUE
)

lenient_by_condition_pair <- input_file_w_filters %>%
    filter(lenient_condition_pair) %>%
    group_by_at(groupby_cols) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = target, values_from = n, values_fill = NA) %>%
    data.frame()
write.xlsx(lenient_by_condition_pair,
    file = number_of_interactions_filename,
    sheetName = "lenient_by_condition_pair", append = TRUE
)

stringent_by_condition <- input_file_w_filters %>%
    filter(stringent_condition) %>%
    group_by_at(groupby_cols) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = target, values_from = n, values_fill = NA) %>%
    data.frame()
write.xlsx(stringent_by_condition,
    file = number_of_interactions_filename,
    sheetName = "stringent_by_condition", append = TRUE
)

stringent_by_condition_pair <- input_file_w_filters %>%
    filter(stringent_condition_pair) %>%
    group_by_at(groupby_cols) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = target, values_from = n, values_fill = NA) %>%
    data.frame()
write.xlsx(stringent_by_condition_pair,
    file = number_of_interactions_filename,
    sheetName = "stringent_by_condition_pair", append = TRUE
)

log_info("Save output...")
saveRDS(input_file_w_filters,
    file = glue("{args$output_dir}/400_samples_interactions_mvoted_w_filters.rds")
)
log_info("COMPLETED!")
