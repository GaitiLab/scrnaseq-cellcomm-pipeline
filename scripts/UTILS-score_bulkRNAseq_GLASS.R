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
    parser$add_argument("--signatures", type = "character", help = "Path to signatures matrix from CIBERSORT")
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
    args$signatures <- "000_misc_local/gene_lists/CIBERSORT_cell_type_signatures.rds"
    args$meta <- "/Users/joankant/Desktop/gaitigroup/Data/GBM/public_data/GLASS/analysis_rnaseq_pairs.csv"
    args$gene_exp <- "/Users/joankant/Desktop/gaitigroup/Data/GBM/public_data/GLASS/gene_tpm_matrix_all_samples.tsv"
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

# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#3_Overview_of_the_GSVA_functionality
log_info("Compute ssGSEA scores for each signature...")
ssgsea <- GSVA::gsva(as(data.matrix(gene_exp_mat), "dgCMatrix"), signatures_list, method = "ssgsea")

log_info("Compute scalop scores for each signature (Neftel et al. 2019)...")
scalop_scores <- sigScores(data.matrix(gene_exp_mat), sigs = signatures_list)

saveRDS(list(ssgsea = ssgsea, scalop = scalop_scores), file = glue("{args$output_dir}/GLASS_paired_scored_{args$suffix}.rds"))
log_info("COMPLETED!")
