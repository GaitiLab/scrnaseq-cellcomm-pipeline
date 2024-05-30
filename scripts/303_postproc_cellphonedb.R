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
    args$sample_id <- "6509_cortex"
    args$output_dir <- glue("{here::here()}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/303_postproc_cpdb")
    args$interaction_scores <- glue("{here::here()}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/203_cci_cpdb/statistical_analysis_interaction_scores__{args$sample_id}.txt")
    args$pval <- glue("{here::here()}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/203_cci_cpdb/statistical_analysis_pvalues__{args$sample_id}.txt")
    args$sign_means <- glue("{here::here()}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/203_cci_cpdb/statistical_analysis_significant_means__{args$sample_id}.txt")
    args$means <- glue("{here::here()}/output/CCI_CellClass_L2_2_reassigned_samples_confident_only/203_cci_cpdb/statistical_analysis_means__{args$sample_id}.txt")
    args$ref_db <- glue("{here::here()}/data/interactions_db/ref_db.rds")
}

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

# Setting sample + run id depending on downsampling or not
split_sample_id <- str_split(args$sample_id, "__", simplify = TRUE)
sample_id <- split_sample_id[1]
run_id <- NA
if (length(split_sample_id) > 2) {
    run_id <- split_sample_id[3]
}

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
    mutate(method = "CellPhoneDBv5", Sample = sample_id, run_id = run_id) %>%
    rename(CellPhoneDB_score = interaction_score)
#   interacting_pair        source_target pval  rank sign_mean      mean interaction_score simple_interaction complex_interaction        method         Sample run_id
# 1         NRG2_MOG Malignant__Malignant    1 1.889        NA 0.7974264                 0          NRG2__MOG           NRG2__MOG CellPhoneDBv5 6234_2895153_A      1
# 2        NRG2_NRP2 Malignant__Malignant    1 0.222        NA 0.7786296                 0         NRG2__NRP2          NRG2__NRP2 CellPhoneDBv5 6234_2895153_A      1
# 3       VEGFD_NRP2 Malignant__Malignant    1 1.889        NA 0.0150624                 0        VEGFD__NRP2         VEGFD__NRP2 CellPhoneDBv5 6234_2895153_A      1
# 4       VEGFA_NRP2 Malignant__Malignant    1 0.222        NA 0.1537503                 0        VEGFA__NRP2         VEGFA__NRP2 CellPhoneDBv5 6234_2895153_A      1
# 5         PGF_NRP2 Malignant__Malignant    1 1.889        NA 0.0295120                 0          PGF__NRP2           PGF__NRP2 CellPhoneDBv5 6234_2895153_A      1
# 6       VEGFC_NRP2 Malignant__Malignant    1 1.889        NA 0.0000000                 0        VEGFC__NRP2         VEGFC__NRP2 CellPhoneDBv5 6234_2895153_A      1

log_info("Save output...")
saveRDS(interactions, glue("{args$output_dir}/cpdb__{args$sample_id}__postproc.rds"))

log_info("COMPLETED!")
