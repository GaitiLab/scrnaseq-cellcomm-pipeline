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
    parser$add_argument("--sample_varname", type = "character", help = "Name of sample variable", default = "")
    parser$add_argument("--condition_varname", type = "character", help = "Name of condition variable", default = "")
    parser$add_argument("--patient_varname", type = "character", help = "Name of patient variable", default = "")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$input_dir <- glue("{here::here()}/output/CCI_CellClass_L1_updated/400_consensus")
    args$output_dir <- glue("{here::here()}/output/CCI_CellClass_L1_updated/401_combine_samples")
    args$celltype_oi <- NULL
    args$metadata <- glue("{here::here()}/output/CCI_CellClass_L1_updated/000_data/gbm_regional_study__metadata.rds")
    args$meta_vars_oi <- glue("{here::here()}/000_misc_local/meta_vars_oi.txt")
    args$sample_varname <- "Sample"
    # args$condition_varname <- "Region_Grouped"
    args$patient_varname  <- "Patient"
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

log_info("Load majority voted interactions...")
all_mvoted <- list.files(args$input_dir,
    full.names = TRUE, pattern = glue("*__interactions_mvoted.rds")
)
all_mvoted <- do.call(rbind, lapply(all_mvoted, readRDS))
# r$> head(all_mvoted)
# # A tibble: 6 × 10
# # Groups:   Sample, source_target, complex_interaction [6]
#   Sample         source_target        complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb lenient_voting stringent_voting
#   <chr>          <chr>                <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl> <lgl>          <lgl>
# 1 6234_2895153_A Malignant__Malignant A2M__LRP1                   1        0           0            1       0 FALSE          FALSE
# 2 6234_2895153_A Malignant__Malignant ACE__BDKRB2                 1        0           0            1       0 FALSE          FALSE
# 3 6234_2895153_A Malignant__Malignant ADAM10__AXL                 1        0           0            1       0 FALSE          FALSE
# 4 6234_2895153_A Malignant__Malignant ADAM10__CADM1               1        0           0            1       0 FALSE          FALSE
# 5 6234_2895153_A Malignant__Malignant ADAM10__GPNMB               1        0           0            1       0 FALSE          FALSE
# 6 6234_2895153_A Malignant__Malignant ADAM10__IL6R                1        0           0            1       0 FALSE          FALSE

# r$> tail(all_mvoted)
# # A tibble: 6 × 10
# # Groups:   Sample, source_target, complex_interaction [6]
#   Sample          source_target        complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb lenient_voting stringent_voting
#   <chr>           <chr>                <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl> <lgl>          <lgl>
# 1 6514_solid_core Microglia__Microglia WNT5A__PTK7                 1        0           0            1       0 FALSE          FALSE
# 2 6514_solid_core Microglia__Microglia WNT5A__PTPRK                1        0           0            1       0 FALSE          FALSE
# 3 6514_solid_core Microglia__Microglia WNT5A__ROR2                 1        0           0            1       0 FALSE          FALSE
# 4 6514_solid_core Microglia__Microglia WNT5A__RYK                  2        0           0            1       1 FALSE          FALSE
# 5 6514_solid_core Microglia__Microglia YBX1__NOTCH1                1        0           0            1       0 FALSE          FALSE
# 6 6514_solid_core Microglia__Microglia ZP3__MERTK                  1        0           0            1       0 FALSE          FALSE

log_info("Load significant interactions (includes sample information)...")
all_sign_interactions <- list.files(args$input_dir,
    full.names = TRUE, pattern = glue("*__signif_interactions.rds")
)
all_sign_interactions <- do.call(rbind, lapply(all_sign_interactions, readRDS))
# r$> head(all_sign_interactions)
# # A tibble: 6 × 8
# # Groups:   Sample, source_target, complex_interaction [6]
#   Sample         source_target        complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb
#   <chr>          <chr>                <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl>
# 1 6234_2895153_A Malignant__Malignant A2M__LRP1                   1        0           0            1       0
# 2 6234_2895153_A Malignant__Malignant ACE__BDKRB2                 1        0           0            1       0
# 3 6234_2895153_A Malignant__Malignant ADAM10__AXL                 1        0           0            1       0
# 4 6234_2895153_A Malignant__Malignant ADAM10__CADM1               1        0           0            1       0
# 5 6234_2895153_A Malignant__Malignant ADAM10__GPNMB               1        0           0            1       0
# 6 6234_2895153_A Malignant__Malignant ADAM10__IL6R                1        0           0            1       0

# r$> tail(all_sign_interactions)
# # A tibble: 6 × 8
# # Groups:   Sample, source_target, complex_interaction [6]
#   Sample          source_target        complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb
#   <chr>           <chr>                <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl>
# 1 6514_solid_core Microglia__Microglia WNT5A__PTK7                 1        0           0            1       0
# 2 6514_solid_core Microglia__Microglia WNT5A__PTPRK                1        0           0            1       0
# 3 6514_solid_core Microglia__Microglia WNT5A__ROR2                 1        0           0            1       0
# 4 6514_solid_core Microglia__Microglia WNT5A__RYK                  2        0           0            1       1
# 5 6514_solid_core Microglia__Microglia YBX1__NOTCH1                1        0           0            1       0
# 6 6514_solid_core Microglia__Microglia ZP3__MERTK                  1        0           0            1       0


log_info("Load metadata...")
if (file.exists(args$meta_vars_oi) && !is.null(args$meta_vars_oi) && args$meta_vars_oi != "") {
    log_info("Metadata variables of interest provided. Using these variables...")
    cols_oi <- read.csv(args$meta_vars_oi,
        header = FALSE, stringsAsFactors = FALSE
    ) %>%
        pull(V1)
} else {
    log_info("No metadata variables of interest provided. Using default variables (Sample, Region_Grouped)...")
    cols_oi <- unique(c(args$sample_varname, args$condition_varname, args$patient_varname))

}
metadata <- readRDS(args$metadata) 
cols <- colnames(metadata)
cols_oi <- intersect(cols_oi, cols)
metadata <- metadata %>%
    select(all_of(cols_oi)) %>%
    rownames_to_column("tmp") %>%
    select(-tmp) %>%
    distinct()
# NOTE: should only contain the information on a sample-level, not the cell-level.
# r$> head(metadata)
#           Sample Region_Grouped
# 1 6237_2222190_A             PT
# 2 6237_2222190_C             TE
# 3 6237_2222190_D             SC
# 4 6237_2222190_F             SC
# 5 6234_2895153_A             TE
# 6 6234_2895153_B             PT

log_info("Add metadata...")
# TODO: maybe make this optional?
all_mvoted <- all_mvoted %>%
    # All post-processed interaction results have the same column names, e.g. "Sample"
    left_join(metadata, by = c("Sample" = args$sample_varname)) %>%
    ungroup()

if (args$sample_varname == args$patient_varname) {
    all_mvoted <- all_mvoted %>% mutate(Patient = Sample)
}

# r$> head(all_mvoted)
# # A tibble: 6 × 11
# # Groups:   Sample, source_target, complex_interaction [6]
#   Sample         source_target        complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb lenient_voting stringent_voting Region_Grouped
#   <chr>          <chr>                <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl> <lgl>          <lgl>            <chr>
# 1 6234_2895153_A Malignant__Malignant A2M__LRP1                   1        0           0            1       0 FALSE          FALSE            TE
# 2 6234_2895153_A Malignant__Malignant ACE__BDKRB2                 1        0           0            1       0 FALSE          FALSE            TE
# 3 6234_2895153_A Malignant__Malignant ADAM10__AXL                 1        0           0            1       0 FALSE          FALSE            TE
# 4 6234_2895153_A Malignant__Malignant ADAM10__CADM1               1        0           0            1       0 FALSE          FALSE            TE
# 5 6234_2895153_A Malignant__Malignant ADAM10__GPNMB               1        0           0            1       0 FALSE          FALSE            TE
# 6 6234_2895153_A Malignant__Malignant ADAM10__IL6R                1        0           0            1       0 FALSE          FALSE            TE
all_sign_interactions <- all_sign_interactions %>%
    # All post-processed interaction results have the same column names, e.g. "Sample"
    left_join(metadata, by = c("Sample" = args$sample_varname)) %>%
    ungroup()
# r$> head(all_sign_interactions)
# # A tibble: 6 × 9
# # Groups:   Sample, source_target, complex_interaction [6]
#   Sample         source_target        complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb Region_Grouped
#   <chr>          <chr>                <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl> <chr>
# 1 6234_2895153_A Malignant__Malignant A2M__LRP1                   1        0           0            1       0 TE
# 2 6234_2895153_A Malignant__Malignant ACE__BDKRB2                 1        0           0            1       0 TE
# 3 6234_2895153_A Malignant__Malignant ADAM10__AXL                 1        0           0            1       0 TE
# 4 6234_2895153_A Malignant__Malignant ADAM10__CADM1               1        0           0            1       0 TE
# 5 6234_2895153_A Malignant__Malignant ADAM10__GPNMB               1        0           0            1       0 TE
# 6 6234_2895153_A Malignant__Malignant ADAM10__IL6R                1        0           0            1       0 TE

log_info("Save...")
saveRDS(all_mvoted, file = glue("{output_dir}/401_samples_interactions_mvoted.rds"))
saveRDS(all_sign_interactions, file = glue("{output_dir}/401_samples_sign_interactions.rds"))
log_info("COMPLETED!")
