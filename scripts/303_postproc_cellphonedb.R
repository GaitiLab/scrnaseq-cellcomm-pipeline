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
        description = "Post-processing CellPhoneDB output",
    )
    parser$add_argument("--sample_id",
        default = "",
        type = "character", help = "Sample ID"
    )
    parser$add_argument("--interaction_scores",
        default = "",
        type = "character", help = "Path to CellPhoneDB interaction scores"
    )
    parser$add_argument("--pval",
        default = "",
        type = "character", help = "Path to CellPhoneDB p-values"
    )
    parser$add_argument("--sign_means",
        default = "",
        type = "character", help = "Path to CellPhoneDB significant means"
    )
    parser$add_argument("--means",
        default = "",
        type = "character", help = "Path to CellPhoneDB means"
    )
    parser$add_argument("--ref_db",
        default = "",
        type = "character", help = "Path to reference database"
    )

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/303_postproc_cpdbv5")
    args$sample_id <- "6234_2895153_A"
    args$interaction_scores <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/203_cci_cpdbv5/statistical_analysis_interaction_scores__{args$sample_id}.txt")
    args$pval <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/203_cci_cpdbv5/statistical_analysis_pvalues__{args$sample_id}.txt")
    args$sign_means <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/203_cci_cpdbv5/statistical_analysis_significant_means__{args$sample_id}.txt")
    args$means <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/203_cci_cpdbv5/statistical_analysis_means__{args$sample_id}.txt")
    args$ref_db <- glue("{here::here()}/001_data_local/interactions_db_v2/ref_db.rds")
}

# output/CellClass_L4_min3_types_rerun/203_cci_cpdbv5

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)


log_info("Load reference database...")
ref_db <- readRDS(args$ref_db) %>% select(interaction, simple_interaction, complex_interaction, interaction)

log_info("Load CellPhoneDB output...")
pval <- read.table(args$pval, sep = "\t", header = TRUE, check.names = FALSE)
sign_means <- read.table(args$sign_means, sep = "\t", header = TRUE, check.names = FALSE)
means <- read.table(args$means, sep = "\t", header = TRUE, check.names = FALSE)
interaction_scores <- read.table(args$interaction_scores, sep = "\t", header = TRUE, check.names = FALSE)

log_info("Process 'pval'...")
colnames_oi <- colnames(pval)[str_detect(colnames(pval), "\\|")]
pval <- pval %>%
    select(interacting_pair, all_of(colnames_oi)) %>%
    reshape2::melt("interacting_pair", value.name = "pval", variable.name = "source_target") %>%
    mutate(source_target = str_replace_all(source_target, "\\|", "__"))

log_info("Process 'sign_means'...")
colnames_oi <- colnames(sign_means)[str_detect(colnames(sign_means), "\\|")]
sign_means <- sign_means %>%
    select(interacting_pair, rank, all_of(colnames_oi)) %>%
    reshape2::melt(c("interacting_pair", "rank"), value.name = "sign_mean", variable.name = "source_target") %>%
    mutate(source_target = str_replace_all(source_target, "\\|", "__"))

log_info("Process 'means'...")
colnames_oi <- colnames(means)[str_detect(colnames(means), "\\|")]
means <- means %>%
    select(interacting_pair, all_of(colnames_oi)) %>%
    reshape2::melt(c("interacting_pair"), value.name = "mean", variable.name = "source_target") %>%
    mutate(source_target = str_replace_all(source_target, "\\|", "__"))

log_info("Process 'interaction_scores'...")
colnames_oi <- colnames(interaction_scores)[str_detect(colnames(interaction_scores), "\\|")]
interaction_scores <- interaction_scores %>%
    select(interacting_pair, all_of(colnames_oi)) %>%
    reshape2::melt(c("interacting_pair"), value.name = "interaction_score", variable.name = "source_target") %>%
    mutate(source_target = str_replace_all(source_target, "\\|", "__"))

log_info("Merge CellPhoneDB output...")
interactions <- pval %>%
    left_join(sign_means, by = c("interacting_pair", "source_target")) %>%
    left_join(means, by = c("interacting_pair", "source_target")) %>%
    left_join(interaction_scores,
        by = c("interacting_pair", "source_target")
    ) %>%
    left_join(ref_db, by = c("interacting_pair" = "interaction")) %>%
    mutate(method = "CellPhoneDBv5", Sample = args$sample_id)

log_info("Save output...")
saveRDS(interactions, glue("{args$output_dir}/cpdb__{args$sample_id}__postproc.rds"))

log_info("COMPLETED!")
