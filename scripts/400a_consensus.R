if (dir.exists("/opt/.renv")) {
    print("running with Docker/Singularity")
    renv::load("/opt/.renv")
}

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
        description = "Determine consensus for a sample",
        default_output = "400_consensus"
    )
    parser$add_argument("-a", "--alpha",
        type = "double",
        default = 0.05, help = "Significance threshold"
    )
    parser$add_argument("-id", "--sample_id",
        type = "character", default = 1,
        help = "Sample id"
    )
    parser$add_argument("--cellchat_obj",
        type = "character", default = NULL,
        help = "CellChat object"
    )
    parser$add_argument("--liana_obj",
        type = "character", default = NULL,
        help = "LIANA object"
    )
    parser$add_argument("--cell2cell_obj",
        type = "character", default = NULL,
        help = "Cell2Cell object"
    )
    parser$add_argument("--cpdb_obj",
        type = "character", default = NULL,
        help = "CPDB object"
    )
    args <- parser$parse_args()
} else {
    run_dir <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun")

    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types_rerun/400_consensus")
    args$Sample <- "6514_enhancing_border"
    args$alpha <- 0.05

    args$cellchat_obj <- glue("{run_dir}/300_postproc_cellchat/{args$Sample}__postproc.rds")
    args$liana_obj <- glue("{run_dir}/301_postproc_liana/{args$Sample}__postproc.rds")
    args$cell2cell_obj <- glue("{run_dir}/302_postproc_cell2cell/{args$Sample}__postproc.rds")
    args$cpdb_obj <- glue("{run_dir}/303_postproc_cpdb/{args$Sample}__postproc.rds")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- paste0(args$output_dir)
create_dir(output_dir)

log_info("Load postprocessed objects from CCIs...")
common_cols <- c("source_target", "complex_interaction", "pval", "method", "Sample")
obj_cellchat <- readRDS(args$cellchat_obj) %>% select(all_of(common_cols))
obj_liana <- readRDS(args$liana_obj) %>% select(all_of(common_cols))
obj_cell2cell <- readRDS(args$cell2cell_obj) %>% select(all_of(common_cols))
obj_cpdb <- readRDS(args$cpdb_obj) %>% select(all_of(common_cols))

log_info(glue("Number of interactions in CellChat BEFORE filtering: {nrow(obj_cellchat)}"))
log_info(glue("Number of interactions in LIANA BEFORE filtering: {nrow(obj_liana)}"))
log_info(glue("Number of interactions in Cell2Cell BEFORE filtering: {nrow(obj_cell2cell)}"))
log_info(glue("Number of interactions in CPDB BEFORE filtering: {nrow(obj_cpdb)}"))

obj_cellchat <- obj_cellchat %>%
    filter(pval < args$alpha) %>%
    select(-pval)
obj_liana <- obj_liana %>%
    filter(pval < args$alpha) %>%
    select(-pval)
obj_cell2cell <- obj_cell2cell %>%
    filter(pval < args$alpha) %>%
    select(-pval)
obj_cpdb <- obj_cpdb %>%
    filter(pval < args$alpha) %>%
    select(-pval)

log_info(glue("Number of interactions in CellChat AFTER filtering: {nrow(obj_cellchat)}"))
log_info(glue("Number of interactions in LIANA AFTER filtering: {nrow(obj_liana)}"))
log_info(glue("Number of interactions in Cell2Cell AFTER filtering: {nrow(obj_cell2cell)}"))
log_info(glue("Number of interactions in CPDB AFTER filtering: {nrow(obj_cpdb)}"))

log_info("Combine the different methods...")
interactions_signif <- do.call(rbind, list(
    obj_cellchat,
    obj_liana,
    obj_cell2cell,
    obj_cpdb
))

log_info(glue("Number of interactions: {nrow(interactions_signif)}"))
utils::capture.output(print(head(interactions_signif)))
# Remove duplicates
interactions_signif <- interactions_signif %>%
    distinct()
log_info(glue("After filtering: {nrow(interactions_signif)}"))

# Number of votes per interaction per source-target per Sample (i.e. detected in how many methods)
cols_oi <- c("Sample", "source_target", "complex_interaction")
all_votes <- interactions_signif %>%
    group_by(Sample, source_target, complex_interaction) %>%
    count() %>%
    rename(n_methods = n) %>%
    left_join(obj_liana %>% select(all_of(cols_oi)) %>% mutate(in_liana = 1)) %>%
    left_join(obj_cellchat %>% select(all_of(cols_oi)) %>% mutate(in_cellchat = 1)) %>%
    left_join(obj_cell2cell %>% select(all_of(cols_oi)) %>% mutate(in_cell2cell = 1)) %>%
    left_join(obj_cpdb %>% select(all_of(cols_oi)) %>% mutate(in_cpdb = 1))
all_votes[is.na(all_votes)] <- 0
head(all_votes)

log_info("Take consensus...")
interactions_mvoted <- all_votes %>%
    group_by(Sample, source_target, complex_interaction) %>%
    mutate(
        lenient_voting = (n_methods >= 2) & (in_liana == 1),
        stringent_voting = (n_methods >= 3) & (in_liana == 1)
    )

log_info(glue("Number of interactions after consensus: {nrow(interactions_mvoted)}"))

log_info("Save output...")
saveRDS(
    interactions_mvoted,
    glue("{output_dir}/{args$sample_id}__interactions_mvoted.rds")
)
saveRDS(
    all_votes,
    glue("{output_dir}/{args$sample_id}__signif_interactions.rds")
)

log_info("COMPLETED!")
