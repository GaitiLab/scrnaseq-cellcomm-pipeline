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
        description = "Post-processing LIANA results",
    )
    parser$add_argument("--input_interactions", type = "character", default = NULL, help = "Path to LIANA results")
    parser$add_argument("--ref_db", type = "character", default = NULL, help = "Path to reference database")
    parser$add_argument("--sample_id", type = "character", default = NULL, help = "Sample ID")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CCI_CellClass_L1_conf_malign/301_postproc_liana")
    args$input_interactions <- glue("{here::here()}/output/CCI_CellClass_L1_conf_malign/201_cci_liana/liana__6234_2895153_A.rds")
    args$ref_db <- glue("{here::here()}/data/interactions_db/ref_db.rds")
    args$sample_id <- "6234_2895153_A"
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
pacman::p_load(reshape2)
pacman::p_load_gh("saezlab/liana")

log_info("Load data...")
interactions <- readRDS(args$input_interactions)

log_info("Load reference database...")
ref_db <- readRDS(args$ref_db) %>%
    select(
        simple_interaction,
        complex_interaction, interaction
    )

# Setting sample + run id depending on downsampling or not
split_sample_id <- str_split(args$sample_id, "__", simplify = TRUE)
sample_id <- split_sample_id[1]
run_id <- NA
if (length(split_sample_id) > 2) {
    run_id <- split_sample_id[3]
}


log_info(glue("Processing sample={args$sample_id}..."))
interactions <- interactions %>%
    liana::liana_aggregate() %>%
    dplyr::rename(
        pval = aggregate_rank
    ) %>%
    mutate(interaction_score = sca.LRscore) %>%
    # select(source, target, ligand.complex, receptor.complex, pval, interaction_score) %>%
    unite(interaction, ligand.complex, receptor.complex, sep = "_") %>%
    left_join(ref_db, by = "interaction") %>%
    select(-interaction) %>%
    unite(source_target, source, target, sep = "__") %>%
    mutate(method = "LIANA", Sample = sample_id, run_id = run_id)
# # A tibble: 6 × 8
#   source_target                 pval interaction_score simple_interaction complex_interaction method Sample         run_id
#   <chr>                        <dbl>             <dbl> <chr>              <chr>               <chr>  <chr>          <chr>
# 1 Microglia__Microglia 0.00000000723             0.901 C3__ITGAX          C3__ITGAX           LIANA  6234_2895153_A 3
# 2 Microglia__Microglia 0.0000000387              0.890 C3__ITGAX          C3__ITGAX:ITGB2     LIANA  6234_2895153_A 3
# 3 Microglia__Microglia 0.0000000433              0.873 C3__ITGB2          C3__ITGB2           LIANA  6234_2895153_A 3
# 4 Microglia__Microglia 0.000000106               0.833 PECAM1__PECAM1     PECAM1__PECAM1      LIANA  6234_2895153_A 3
# 5 Microglia__Microglia 0.000000121               0.893 C3__NRP1           C3__NRP1            LIANA  6234_2895153_A 3
# 6 Microglia__Microglia 0.000000137               0.785 ICAM1__ITGAX       ICAM1__ITGAX        LIANA  6234_2895153_A 3

log_info("Save output...")
saveRDS(interactions, glue("{args$output_dir}/liana__{args$sample_id}__postproc.rds"))
log_info("COMPLETED!")
