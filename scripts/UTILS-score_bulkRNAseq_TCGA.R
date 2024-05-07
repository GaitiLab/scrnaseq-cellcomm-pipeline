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
        description = "Compute scores for bulkRNAseq TCGA",
    )
    parser$add_argument("--signatures", type = "character", help = "Path to signatures matrix from CIBERSORT")
    parser$add_argument("--meta", type = "character", help = "Path to meta file")
    parser$add_argument("--gene_exp", type = "character", help = "File to gene expression matrix")
    parser$add_argument("--status_info", type = "character", help = "IDH status information")
    parser$add_argument("--suffix", type = "character", help = "Suffix to add for filename", default = "")

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/corr_analysis_bulk_rnaseq")
    args$signatures <- "000_misc_local/gene_lists/CIBERSORT_cell_type_signatures.rds"
    args$interactions <- "/Users/joankant/Library/CloudStorage/OneDrive-UHN/Spatial_GBM/Analysis/WIP/scRNAseq/CCI/unique_interactions.xlsx"
    args$meta <- "/Users/joankant/Desktop/gaitigroup/Data/GBM/public_data/TCGA_bulkRNAseq/clinical.project-tcga-gbm.2024-04-17/clinical.tsv"
    args$gene_exp <- "/Users/joankant/Desktop/gaitigroup/Data/GBM/public_data/TCGA_bulkRNAseq/gdac.broadinstitute.org_GBM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/GBM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt"
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
pacman::p_load(GSVA, boot, scalop, readxl)

log_info("Load gene-expression...")
gene_exp <- read.table(args$gene_exp, check.names = FALSE, sep = "\t", header = TRUE, row.names = 1)

gene_exp_mat <- gene_exp[, gene_exp["gene_id", ] == "scaled_estimate"]
gene_exp_mat <- gene_exp_mat[-1, ]

# Extract gene symbols from Hybridization Ref (rownames)
genes <- str_split(rownames(gene_exp_mat), "\\|", simplify = TRUE)[, 1]

# Remove rows without gene symbols
gene_exp_mat <- gene_exp_mat[genes != "?", ]

# Convert all columns to numeric data type
gene_exp_mat <- gene_exp_mat %>%
    mutate_if(is.character, as.numeric) %>%
    mutate(gene = str_split(rownames(gene_exp_mat), "\\|", simplify = TRUE)[, 1]) %>%
    remove_rownames() %>%
    group_by(gene) %>%
    summarize(across(where(is.numeric), mean))

gene_exp_mat <- gene_exp_mat %>% column_to_rownames("gene")

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

# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#3_Overview_of_the_GSVA_functionality
log_info("Compute ssGSEA scores for each signature...")
ssgsea <- GSVA::gsva(as(data.matrix(gene_exp_mat), "dgCMatrix"), signatures_list, method = "ssgsea")


log_info("Compute scalop scores for each signature (Neftel et al. 2019)...")
scalop_scores <- sigScores(data.matrix(gene_exp_mat), sigs = signatures_list)

# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#3_Overview_of_the_GSVA_functionality
log_info("Compute ssGSEA scores for each signature...")
ssgsea <- GSVA::gsva(as(data.matrix(gene_exp_mat), "dgCMatrix"), signatures_list, method = "ssgsea")

log_info("Compute scalop scores for each signature (Neftel et al. 2019)...")
scalop_scores <- sigScores(data.matrix(gene_exp_mat), sigs = signatures_list)

saveRDS(list(ssgsea = ssgsea, scalop = scalop_scores), file = glue("{args$output_dir}/TCGA_scored_{args$suffix}.rds"))
log_info("COMPLETED!")


log_info("COMPLETED!")
