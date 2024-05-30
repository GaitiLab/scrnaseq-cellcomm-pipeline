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
        description = "Aggregate interactions by sample (combine p-values)",
    )
    parser$add_argument("--input_file", type = "character", help = "Input from 401_combine_samples.R, name should be 401_samples_interactions_agg_rank.rds", default = "")
    parser$add_argument("--condition_var", default = "Region_Grouped", type = "character", help = "Grouping variable")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/GBM_CCI_Analysis/output/CCI_CellClass_L2_2_reassigned_samples_confident_only_FINAL/402_aggregation")
    args$input_file <- glue("/Users/joankant/Desktop/gaitigroup/Users/Joan/GBM_CCI_Analysis/output/CCI_CellClass_L2_2_reassigned_samples_confident_only_FINAL/401_combine_samples/401_samples_interactions_agg_rank.rds")
    args$condition_var <- "Region"
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
log_info("Load additional libraries...")
pacman::p_load(survcomp)

log_info("Load interactions and combine p-values...")
obj <- readRDS(args$input_file) %>%
    ungroup() %>%
    # TODO test, If score missing assign 0 (higher scores = higher interaction scores)
    mutate(
        LIANA_score = ifelse(is.na(LIANA_score), 0, LIANA_score),
        CellPhoneDB_score = ifelse(is.na(CellPhoneDB_score), 0, CellPhoneDB_score),
        CellChat_score = ifelse(is.na(CellChat_score), 0, CellChat_score)
    ) %>%
    select(!!sym(args$condition_var), complex_interaction, pval, LIANA_score, CellPhoneDB_score, CellChat_score, source_target) %>%
    group_by(!!sym(args$condition_var), source_target, complex_interaction) %>%
    # Combine p-values and interaction scores (CellChat, LIANA, CellphoneDB) across samples
    mutate(
        pval = combine.test(pval),
        LIANA_score = mean(LIANA_score),
        CellPhoneDB_score = mean(CellPhoneDB_score),
        CellChat_score = mean(CellChat_score)
    )
# r$> head(obj)
# # A tibble: 6 × 7
# # Groups:   Region, source_target, complex_interaction [6]
#   Region complex_interaction          pval LIANA_score CellPhoneDB_score CellChat_score source_target
#   <chr>  <chr>                       <dbl>       <dbl>             <dbl>          <dbl> <chr>
# 1 PT     PTN__PTPRZ1         0.00000000490       0.918              90.9       0.00295  Invasive-high OPC/NPC1__Invasive-high OPC/NPC1
# 2 PT     BCAN__NRCAM         0.0000000126        0.895             100         0.00214  Invasive-high OPC/NPC1__Invasive-high OPC/NPC1
# 3 PT     NCAM1__PTPRZ1       0.0000000203        0.929             100         0.00340  Invasive-high OPC/NPC1__Invasive-high OPC/NPC1
# 4 PT     SEMA5A__PLXNB3      0.0000000507        0.838             100         0.000458 Invasive-high OPC/NPC1__Invasive-high OPC/NPC1
# 5 PT     CNTN1__PTPRZ1       0.0000000554        0.931             100         0.00353  Invasive-high OPC/NPC1__Invasive-high OPC/NPC1
# 6 PT     DSCAM__DCC          0.0000000760        0.917             100         0.00296  Invasive-high OPC/NPC1__Invasive-high OPC/NPC1

log_info("Save results...")
saveRDS(obj, glue("{args$output_dir}/402b_aggregation_samples.rds"))

log_info("COMPLETED!")
