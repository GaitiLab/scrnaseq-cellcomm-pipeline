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
        description = "Create all figures for GBM project",
    )
    parser$add_argument("--interactions_run_dir", type = "character", help = "Directory with interactions")
    parser$add_argument("--condition_varname", type = "character", help = "Directory with interactions")
    parser$add_argument("--unique_interactions", type = "character", help = "Excel file with unique interactions")
    parser$add_argument("--top_n", type = "numeric", help = "Top-N to label/plot", default = 10)
    parser$add_argument("--source_oi", type = "character", help = "Source cell type of interest", default = "")
    parser$add_argument("--target_oi", type = "character", help = "Target cell type of interest", default = "")
    parser$add_argument("--condition_oi", type = "character", help = "Condition of interest, e.g. region, mutation")
    parser$add_argument("--interactions_db", type = "character", help = "Database with interactions (reference)")
    parser$add_argument("--meta", type = "character", help = "Metadata file (RDS)")
    parser$add_argument("--avg_expr", type = "character", help = "Average expression by Sample and cell type label")

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$run_name <- "CCI_CellClass_L2_2_reassigned_samples_confident_only"
    args$output_dir <- "/Users/joankant/Library/CloudStorage/OneDrive-UHN/Spatial_GBM/Analysis/WIP/scRNAseq/CCI/Figures"
    args$interactions_run_dir <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/scrnaseq-cellcomm/output/{args$run_name}"
    args$condition_varname <- "Region"
    args$unique_interactions <- "/Users/joankant/Library/CloudStorage/OneDrive-UHN/Spatial_GBM/Analysis/WIP/scRNAseq/CCI/{args$run_name}_unique_interactions_neuron_invasive_high.xlsx"
    args$top_n <- 10
    args$source_oi <- "Neuron"
    args$target_oi <- "Invasive-high OPC/NPC1"
    args$condition_oi <- "PT"
    args$interactions_db <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/scrnaseq-cellcomm/data/interactions_db/ref_db.rds"
    args$meta <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/scrnaseq-cellcomm/output/{args$run_name}/000_data/gbm_regional_study__metadata.rds"
    args$avg_expr <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/scrnaseq-cellcomm/output/average_expression_and_pseudobulk/average_expression_by_Sample_CCI_CellClass_L2_2.rds"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Manual
args$source_oi <- "Neuron"
args$target_oi <- "Invasive-high OPC/NPC1"

# Main figures
log_info("Differential Heatmap: PT and TC...")
rmarkdown::render(glue("{here::here()}/notebooks/FIGURES-num_interactions_diff_undirected.rmd"),
    params = list(
        input_file = paste0(args$interactions_run_dir, "/402_aggregation/402c_aggregation_integration.rds"),
        output_dir = paste0(args$output_dir, "/", get_name(args$interactions_run_dir)),
        condition_varname = args$condition_varname,
        remove_autocrine = 1,
        group1 = "PT",
        group2 = "TC"
    )
)

log_info("Ranking interactions (hockeyplot)...")
rmarkdown::render(glue("{here::here()}/notebooks/FIGURES-ranking_interactions.rmd"),
    params = list(
        input_file = paste0(args$interactions_run_dir, "/402_aggregation/402c_aggregation_integration.rds"),
        unique_interactions = args$unique_interactions,
        output_dir = paste0(args$output_dir, "/", get_name(args$interactions_run_dir)),
        condition_varname = args$condition_varname,
        top_n = args$top_n,
        source = args$source_oi,
        target = args$target_oi,
        condition_oi = args$condition_oi
    )
)

log_info("Heatmap differences between two cell type pairs...")
rmarkdown::render(glue("{here::here()}/notebooks/FIGURES-interactions_diff_cell_type_pairs.rmd"),
    params = list(
        input_file = paste0(args$interactions_run_dir, "/402_aggregation/402c_aggregation_integration.rds"),
        unique_interactions = args$unique_interactions,
        output_dir = paste0(args$output_dir, "/", get_name(args$interactions_run_dir)),
        score = "CellPhoneDB_score",
        top_n = 20
    )
)


log_info("Gene expression scatterplot for interactions of interest...")
rmarkdown::render(glue("{here::here()}/notebooks/FIGURES-LR_avg_gene_expr.rmd"),
    params = list(
        input_file = paste0(args$interactions_run_dir, "/402_aggregation/402c_aggregation_integration.rds"),
        unique_interactions = args$unique_interactions,
        output_dir = paste0(args$output_dir, "/", get_name(args$interactions_run_dir)),
        condition_varname = args$condition_varname,
        top_n = args$top_n,
        source = args$source_oi,
        target = args$target_oi,
        condition_oi = args$condition_oi,
        interactions_db = args$interactions_db,
        cutoff = 3,
        meta = args$meta,
        avg_expr = args$avg_expr
    )
)

# ---- Supplementary ---- #
log_info("SUPP: Heatmap w/ number of interactions (individual)...")
rmarkdown::render(glue("{here::here()}/notebooks/FIGURES-num_interactions_indiv_undirected.rmd"),
    params = list(
        input_file = paste0(args$interactions_run_dir, "/402_aggregation/402c_aggregation_integration.rds"),
        output_dir = paste0(args$output_dir, "/", get_name(args$interactions_run_dir)),
        condition_varname = args$condition_varname,
        remove_autocrine = 1
    )
)

log_info("SUPP: Analysis impact of cells...")
rmarkdown::render(glue("{here::here()}/notebooks/SUPP-analysis_impact_of_cells.rmd"),
    params = list(
        meta = args$meta,
        input_file = paste0(args$interactions_run_dir, "/401_combine_samples/401_samples_interactions_mvoted.rds"),
        output_dir = paste0(args$output_dir, "/", get_name(args$interactions_run_dir)),
        percentile = 0.9,
        annot = "CCI_CellClass_L2",
        condition_varname = args$condition_varname,
        avg_expr = args$avg_expr
    )
)


# Dotplot
log_info("SUPP: Dotplot")
rmarkdown::render(glue("{here::here()}/notebooks/FIGURES-num_interactions_all_regions_undirected.rmd"),
    params = list(
        input_file = paste0(args$interactions_run_dir, "/402_aggregation/402c_aggregation_integration.rds"),
        output_dir = paste0(args$output_dir, "/", get_name(args$interactions_run_dir)),
        condition_varname = args$condition_varname
    )
)


log_info("COMPLETED!")
