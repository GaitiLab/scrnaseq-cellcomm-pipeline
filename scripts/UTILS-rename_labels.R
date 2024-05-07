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
        description = "Rename labels",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
    args$run_name <- "CCI_CellClass_L2_2"
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
combi_agg <- readRDS(glue("output/{args$run_name}/402_aggregation/402c_aggregation_integration.rds"))
#  [1] "Region_Grouped"                 "complex_interaction"            "source_target"                  "pval"                           "lenient_condition_patients"
#  [6] "lenient_condition_n_patients"   "lenient_condition_n_samples"    "lenient_condition_samples"      "lenient_condition"              "stringent_condition_patients"
# [11] "stringent_condition_n_patients" "stringent_condition_n_samples"  "stringent_condition_samples"    "stringent_condition"

new_label <- c("Neuronal OPC-like" = "Invasive-high OPC/NPC1")
new_region_label <- c("SC" = "TC")

combi_agg <- combi_agg %>%
    mutate(source_target = str_replace_all(source_target, pattern = new_label), Region_Grouped = str_replace_all(Region_Grouped, pattern = new_region_label)) %>%
    mutate(Region_Grouped = factor(Region_Grouped, levels = names(GBMutils::load_color_palette("Region"))))


saveRDS(combi_agg, glue("output/{args$run_name}/402_aggregation/402c_aggregation_integration.rds"))
