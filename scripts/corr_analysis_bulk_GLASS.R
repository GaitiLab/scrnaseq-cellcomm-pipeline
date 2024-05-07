# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
require(GBMutils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Bulk RNAseq data Correlation Analysis for GLASS Cohort",
    )
    parser$add_argument("--n_iter", type = "numeric", help = "Number of iterations for bootstrapping", default = 100)
    parser$add_argument("--corr_method", type = "character", help = "Method to use for correlation, e.g. pearson or spearman", default = "spearman")
    parser$add_argument("--signatures", type = "character", help = "Path to signatures matrix from CIBERSORT")
    parser$add_argument("--n_genes", type = "numeric", help = "Number of top-genes to use", default = 5)
    parser$add_argument("--interactions", type = "character", help = "Excel sheet with interactions")
    parser$add_argument("--meta", type = "character", help = "Path to meta file")
    parser$add_argument("--gene_exp", type = "character", help = "File to gene expression matrix")
    parser$add_argument("--status_info", type = "character", help = "IDH status information")
    parser$add_argument("--suffix", type = "character", help = "Suffix to add to output filename", default = "")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/corr_analysis_bulk_rnaseq")
    args$n_iter <- 100
    args$corr_method <- "spearman"
    args$signatures <- "000_misc_local/gene_lists/CIBERSORT_cell_type_signatures.rds"
    args$n_genes <- 5
    args$interactions <- "/Users/joankant/Library/CloudStorage/OneDrive-UHN/Spatial_GBM/Analysis/WIP/scRNAseq/CCI/unique_interactions.xlsx"
    args$meta <- "/Users/joankant/Desktop/gaitigroup/Data/GBM/public_data/GLASS/analysis_rnaseq_pairs.csv"
    args$gene_exp <- "/Users/joankant/Desktop/gaitigroup/Data/GBM/public_data/GLASS/gene_tpm_matrix_all_samples.tsv"
    args$status_info <- "/Users/joankant/Desktop/gaitigroup/Data/GBM/public_data/GLASS/clinical_surgerie.csv"
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
pacman::p_load(GSVA, readxl, boot, scalop)

log_info("Load gene-expression...")
gene_exp_mat <- read.table(args$gene_exp, header = TRUE, row.names = 1, check.names = FALSE)

log_info("Load metadata...")
samples_wt <- read.csv(args$status_info) %>%
    filter(idh_codel_subtype == "IDHwt", !is.na(sample_barcode)) %>%
    pull(case_barcode)

meta <- read.csv(args$meta) %>% filter(case_barcode %in% samples_wt, sample_type_a == "TP", sample_type_b == "R1")

log_info("Load cell type gene signatures...")
signatures_list <- readRDS(args$signatures)
if (str_detect(args$signatures, "CIBERSORT")) {
    # TODO Change if necessary
    celltypes_oi <- c("Malignant_AC", "Malignant_MES_AST", "Malignant_MES_HYP", "Malignant_MES_INT", "Malignant_NPC1", "Malignant_NPC2", "Malignant_OPC", "Neuronal.OPC.like", "Neuron")
    signatures_list <- signatures_list[celltypes_oi]

    # Collapse signatures
    mes_like <- c("Malignant_MES_AST", "Malignant_MES_HYP", "Malignant_MES_INT")
    npc_like <- c("Malignant_NPC1", "Malignant_NPC2")

    mes_like_signature <- signatures_list[mes_like] %>%
        unlist() %>%
        unique()
    npc_like_signature <- signatures_list[npc_like] %>%
        unlist() %>%
        unique()

    # Remove subtypes for mes and npc
    signatures_list[c(mes_like, npc_like)] <- NULL

    signatures_list[["MES.like"]] <- mes_like_signature
    signatures_list[["NPC.like"]] <- npc_like_signature

    # Change naming of OPC and NPC
    names(signatures_list)[names(signatures_list) == "Malignant_AC"] <- "AC.like"
    names(signatures_list)[names(signatures_list) == "Malignant_OPC"] <- "OPC.like"
}

log_info("Selecting samples for Primary vs Recurrent samples...")
glass_primary <- meta %>%
    select(tumor_barcode_a) %>%
    mutate(group = "Primary")


glass_recurrent <- meta %>%
    select(tumor_barcode_b) %>%
    mutate(group = "Recurrent")
colnames(glass_primary) <- c("Sample", "Group")
colnames(glass_recurrent) <- c("Sample", "Group")

glass_paired <- rbind(glass_primary, glass_recurrent)
glass_paired %>%
    pull(Group) %>%
    table()

log_info("Load interactions of interest (order by pval)...")
interactions_subset <- read_excel(args$interactions) %>% arrange(pval)
head(interactions_subset)

genes_oi <- interactions_subset %>%
    separate(complex_interaction, into = c("ligand_complex", "receptor_complex"), sep = "__") %>%
    pull(ligand_complex, receptor_complex) %>%
    str_split(., ":") %>%
    unlist() %>%
    unique()

log_info(glue("Number of unique genes in interactions: {length(genes_oi)}..."))

# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#3_Overview_of_the_GSVA_functionality
log_info("Compute ssGSEA scores for each signature...")
ssgsea <- GSVA::gsva(as(data.matrix(gene_exp_mat), "dgCMatrix"), signatures_list, method = "ssgsea")

log_info("Compute scalop scores for each signature (Neftel et al. 2019)...")
scalop_scores <- sigScores(data.matrix(gene_exp_mat), sigs = signatures_list)

log_info("Subset expression matrix for receptor genes of interest...")
common_genes <- intersect(genes_oi, rownames(gene_exp_mat))
log_info(glue("Number of genes in common: {length(common_genes)}..."))

genes_oi_exp <- gene_exp_mat[common_genes, ]

mat <- cbind(t(data.matrix(genes_oi_exp)), t(data.matrix(ssgsea)))
rownames(mat) <- str_replace_all(rownames(mat), "\\.", "-")
df <- data.frame(mat) %>% rownames_to_column("Sample")
glass_paired_df <- glass_paired %>% left_join(df, by = "Sample")

mat <- cbind(t(data.matrix(genes_oi_exp)), scalop_scores)
rownames(mat) <- str_replace_all(rownames(mat), "\\.", "-")
df <- data.frame(mat) %>% rownames_to_column("Sample")
glass_paired_scalop_df <- glass_paired %>% left_join(df, by = "Sample")

samples_primary_ssgsea <- glass_paired_df %>% filter(Group == "Primary")
samples_recurrent_ssgsea <- glass_paired_df %>% filter(Group == "Recurrent")

samples_primary_scalop <- glass_paired_scalop_df %>% filter(Group == "Primary")
samples_recurrent_scalop <- glass_paired_scalop_df %>% filter(Group == "Recurrent")

# Pairs to test, all malignant state x receptor of interest combinations.
pairs_to_test <- expand.grid(names(signatures_list), common_genes, stringsAsFactors = FALSE)

log_info("Correlation analysis using bootstrapping...")
samples_primary_boot_ssgsea <- do.call(
    rbind,
    apply(pairs_to_test, 1, run_bootstrap_corr, df = samples_primary_ssgsea, stat_func = compute_corr, n_iter = args$n_iter, method = args$corr_method)
) %>% mutate(Group = "Primary")

samples_recurrent_boot_ssgsea <- do.call(
    rbind,
    apply(pairs_to_test, 1, run_bootstrap_corr, df = samples_recurrent_ssgsea, stat_func = compute_corr, n_iter = args$n_iter, method = args$corr_method)
) %>% mutate(Group = "Recurrent")


samples_primary_boot_scalop <- do.call(
    rbind,
    apply(pairs_to_test, 1, run_bootstrap_corr, df = samples_primary_scalop, stat_func = compute_corr, n_iter = args$n_iter, method = args$corr_method)
) %>% mutate(Group = "Primary")

samples_recurrent_boot_scalop <- do.call(
    rbind,
    apply(pairs_to_test, 1, run_bootstrap_corr, df = samples_recurrent_scalop, stat_func = compute_corr, n_iter = args$n_iter, method = args$corr_method)
) %>% mutate(Group = "Recurrent")

# Combine primary and recurrent results
glass_boot_ssgsea <- rbind(samples_primary_boot_ssgsea, samples_recurrent_boot_ssgsea)
glass_boot_scalop <- rbind(samples_primary_boot_scalop, samples_recurrent_boot_scalop)

saveRDS(list(ssgsea = glass_boot_ssgsea, scalop = glass_boot_scalop), file = glue("{args$output_dir}/GLASS_paired_corr_boot_{args$suffix}.rds"))

log_info("COMPLETED!")
