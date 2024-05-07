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
        description = "Extract cell type signatures from CIBERSORT's signature matrix",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/000_misc_local/")
    args$signatures <- "/Users/joankant/Desktop/gaitigroup/Users/Pengcheng/workspace/scGBM/cibersort/20240429/CIBERSORTx_sigmatrix_Adjusted.txt"
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
log_info("Load cell type gene signatures...")
signatures <- read.csv(args$signatures, sep = "\t") %>%
    # filter(GeneSymbol %in% rownames(gene_exp_mat)) %>%
    column_to_rownames(("GeneSymbol"))

log_info("Obtain gene signatures for each cell type without filtering (value > 0)...")
signatures_list <- apply(signatures, 2, function(col) {
    rownames(signatures)[col > 0]
})
log_info("Save in a single object...")
saveRDS(signatures_list, file = glue("{args$output_dir}/CIBERSORT_cell_type_signatures_20240429.rds"))




# # Create CellClass L2
# malignant_mes <- c("Malignant_MES_AST", "Malignant_MES_HYP", "Malignant_MES_INT")
# malignant_npc <- c("Malignant_NPC1", "Malignant_NPC2")

# signatures_mes <- signatures_list[malignant_mes] %>%
#     unlist() %>%
#     unique()
# signatures_npc <- signatures_list[malignant_npc] %>%
#     unlist() %>%
#     unique()

# signatures_list[malignant_mes] <- NULL
# signatures_list[malignant_npc] <- NULL

# to_add <- list(Malignant_MES = signatures_mes, Malignant_NPC = signatures_npc)

# signatures_list <- c(signatures_list, to_add)


# log_info("Save in a single object...")
# saveRDS(signatures_list_filtered, file = glue("{args$output_dir}/CIBERSORT_cell_type_signatures_filtered.rds"))
# saveRDS(signatures_list, file = glue("{args$output_dir}/CIBERSORT_cell_type_signatures_L2.rds"))


log_info("COMPLETED!")
