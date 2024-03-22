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
        type = "numeric",
        default = 0.05, help = "Significance threshold"
    )
    parser$add_argument("-id", "--sample_id",
        type = "character", default = 1,
        help = "Sample id"
    )
    parser$add_argument("--cellchat_obj",
        type = "character", default = "",
        help = "CellChat object"
    )
    parser$add_argument("--liana_obj",
        type = "character", default = "",
        help = "LIANA object"
    )
    parser$add_argument("--cell2cell_obj",
        type = "character", default = "",
        help = "Cell2Cell object"
    )
    parser$add_argument("--cpdb_obj",
        type = "character", default = "",
        help = "CPDB object"
    )
    parser$add_argument("--run_dir",
        type = "character", default = "",
        help = "Path to run directory"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$run_dir <- glue("{here::here()}/output/CCI_CellClass_L1_conf_malign")

    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CCI_CellClass_L1_conf_malign/400_consensus")
    args$sample_id <- "6234_2895153_A"
    args$alpha <- 0.05
    # args$cellchat_obj <- glue("{run_dir}/300_postproc_cellchat/cellchat__{args$sample_id}__postproc.rds")
    # args$liana_obj <- glue("{run_dir}/301_postproc_liana/liana__{args$sample_id}__postproc.rds")
    # args$cell2cell_obj <- glue("{run_dir}/302_postproc_cell2cell/cell2cell__{args$sample_id}__postproc.rds")
    # args$cpdb_obj <- glue("{run_dir}/303_postproc_cpdb/cpdb__{args$sample_id}__postproc.rds")
    args$cellchat_obj <- ""
    args$liana_obj <- ""
    args$cell2cell_obj <- ""
    args$cpdb_obj <- ""
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

pacman::p_load(RobustRankAggreg)


log_info("Load postprocessed objects from CCIs...")
common_cols <- c("source_target", "complex_interaction", "pval", "method", "Sample")
# Pre-filtering on p-value, only impacts when downsampling is done (pval == NA, if interaction not detected in each downsampling run)
obj_cellchat <- readRDS(ifelse(file.exists(args$cellchat_obj),
    args$cellchat_obj,
    glue("{args$run_dir}/300_postproc_cellchat/cellchat__{args$sample_id}__postproc.rds")
))
# select(all_of(common_cols)) %>%
# filter(!is.na(pval))
obj_liana <- readRDS(ifelse(file.exists(args$liana_obj),
    args$liana_obj,
    glue("{args$run_dir}/301_postproc_liana/liana__{args$sample_id}__postproc.rds")
))
# select(all_of(common_cols)) %>%
# filter(!is.na(pval))
obj_cell2cell <- readRDS(ifelse(file.exists(args$cell2cell_obj), args$cell2cell_obj, glue("{args$run_dir}/302_postproc_cell2cell/cell2cell__{args$sample_id}__postproc.rds")))
# select(all_of(common_cols)) %>%
# filter(!is.na(pval))
obj_cpdb <- readRDS(ifelse(file.exists(args$cpdb_obj), args$cpdb_obj, glue("{args$run_dir}/303_postproc_cpdb/cpdb__{args$sample_id}__postproc.rds")))
# select(all_of(common_cols)) %>%
# filter(!is.na(pval))

log_info(glue("Number of interactions in CellChat BEFORE filtering: {nrow(obj_cellchat)}"))
log_info(glue("Number of interactions in LIANA BEFORE filtering: {nrow(obj_liana)}"))
log_info(glue("Number of interactions in Cell2Cell BEFORE filtering: {nrow(obj_cell2cell)}"))
log_info(glue("Number of interactions in CPDB BEFORE filtering: {nrow(obj_cpdb)}"))


# ---- Order from most-important (significant) to least based on interaction
# Sorting by score/specificity (NOT p-value)
connectome <- obj_liana %>%
    select(source_target, complex_interaction, Sample, connectome.weight_sc) %>%
    mutate(method = "Connectome (LIANA)") %>%
    rename(score = connectome.weight_sc) %>%
    group_by(Sample, source_target, ) %>%
    arrange(Sample, source_target, desc(score), .by_group = TRUE)

logfc <- obj_liana %>%
    select(source_target, complex_interaction, Sample, logfc.logfc_comb) %>%
    mutate(method = "LogFC (LIANA)") %>%
    rename(score = logfc.logfc_comb) %>%
    group_by(Sample, source_target, ) %>%
    arrange(Sample, source_target, desc(score), .by_group = TRUE)

natmi <- obj_liana %>%
    select(source_target, complex_interaction, Sample, natmi.edge_specificity) %>%
    mutate(method = "NATMI (LIANA)") %>%
    rename(score = natmi.edge_specificity) %>%
    group_by(Sample, source_target, ) %>%
    arrange(Sample, source_target, desc(score), .by_group = TRUE)

sca <- obj_liana %>%
    select(source_target, complex_interaction, Sample, sca.LRscore) %>%
    mutate(method = "SCA (LIANA)") %>%
    rename(score = sca.LRscore) %>%
    group_by(Sample, source_target, ) %>%
    arrange(Sample, source_target, desc(score), .by_group = TRUE)

cytotalk <- obj_liana %>%
    select(source_target, complex_interaction, Sample, run_id, cytotalk.crosstalk_score) %>%
    mutate(method = "Cytotalk (LIANA)") %>%
    rename(score = cytotalk.crosstalk_score) %>%
    group_by(Sample, source_target, ) %>%
    arrange(Sample, source_target, desc(score), .by_group = TRUE)

# Sorting by p-values
obj_cellchat <- obj_cellchat %>%
    group_by(Sample, source_target, ) %>%
    arrange(Sample, source_target, pval, .by_group = TRUE)

obj_cell2cell <- obj_cell2cell %>%
    group_by(Sample, source_target, ) %>%
    arrange(Sample, source_target, pval, .by_group = TRUE)

obj_cpdb <- obj_cpdb %>%
    group_by(Sample, source_target, ) %>%
    arrange(Sample, source_target, pval, .by_group = TRUE)

# Grab all source-target pairs
source_target_pairs <- unique(c(obj_cellchat %>% pull(source_target) %>% unique(), obj_liana %>% pull(source_target) %>% unique(), obj_cell2cell %>% pull(source_target) %>% unique(), obj_cpdb %>% pull(source_target) %>% unique()))


aggregated_ranks <- lapply(source_target_pairs, rank_interactions, connectome = connectome, sca = sca, cytotalk = cytotalk, natmi = natmi, logfc = logfc, obj_cellchat = obj_cellchat, obj_cell2cell = obj_cell2cell, obj_cpdb = obj_cpdb)

interactions_df_aggregated <- do.call(rbind, aggregated_ranks) %>%
    rename(complex_interaction = Name, pval = Score) %>%
    mutate(Sample = obj_liana %>% pull(Sample) %>% unique()) %>%
    remove_rownames() %>%
    # Adding individual p-values
    left_join(obj_cellchat %>% rename(cellchat = pval) %>% select(Sample, source_target, complex_interaction, cellchat)) %>%
    left_join(obj_cell2cell %>% rename(cell2cell = pval) %>% select(Sample, source_target, complex_interaction, cell2cell)) %>%
    left_join(obj_cell2cell %>% rename(liana = pval) %>% select(Sample, source_target, complex_interaction, liana)) %>%
    left_join(obj_cpdb %>% rename(cpdb = pval) %>% select(Sample, source_target, complex_interaction, cpdb))

#  ===== APPROACH 2: FILTERING BASED ON SIGNIFICANCE ===== #
log_info(glue("Filter by significance... (alpha = {args$alpha})"))
obj_cellchat <- obj_cellchat %>%
    select(all_of(common_cols)) %>%
    filter(!is.na(pval)) %>%
    filter(pval < args$alpha) %>%
    select(-pval)
# r$> head(obj_cellchat)
#                 source_target complex_interaction     method                Sample
# 1        Malignant__Malignant COL4A2__ITGAV:ITGB8 CellChatv2 6514_enhancing_border
# 2  Oligodendrocyte__Malignant COL4A5__ITGAV:ITGB8 CellChatv2 6514_enhancing_border
# 3      Macrophage__Macrophage   SPP1__ITGAV:ITGB1 CellChatv2 6514_enhancing_border
# 4       Microglia__Macrophage   SPP1__ITGAV:ITGB1 CellChatv2 6514_enhancing_border
# 5 Oligodendrocyte__Macrophage   SPP1__ITGAV:ITGB1 CellChatv2 6514_enhancing_border
# 6        Malignant__Malignant   GLS:SLC1A2__GRIA2 CellChatv2 6514_enhancing_border
obj_liana <- obj_liana %>%
    select(all_of(common_cols)) %>%
    filter(!is.na(pval)) %>%
    filter(pval < args$alpha) %>%
    select(-pval)
# r$> head(obj_liana)
# # A tibble: 6 × 4
#   source_target                    complex_interaction method Sample
#   <chr>                            <chr>               <chr>  <chr>
# 1 Oligodendrocyte__Oligodendrocyte CNTN2__CNTN2        LIANA  6514_enhancing_border
# 2 Malignant__Malignant             BCAN__EGFR          LIANA  6514_enhancing_border
# 3 Oligodendrocyte__Oligodendrocyte CLDN11__CLDN11      LIANA  6514_enhancing_border
# 4 Malignant__Malignant             BCAN__NRCAM         LIANA  6514_enhancing_border
# 5 Oligodendrocyte__Oligodendrocyte TF__LRP2            LIANA  6514_enhancing_border
# 6 Malignant__Malignant             CHL1__CHL1          LIANA  6514_enhancing_border
obj_cell2cell <- obj_cell2cell %>%
    select(all_of(common_cols)) %>%
    filter(!is.na(pval)) %>%
    filter(pval < args$alpha) %>%
    select(-pval)
# r$> head(obj_cell2cell)
#           source_target  complex_interaction    method                Sample
# 1 Macrophage__Malignant    NRG2__ERBB2:ERBB3 Cell2Cell 6514_enhancing_border
# 2 Macrophage__Malignant    NRG2__ERBB2:ERBB4 Cell2Cell 6514_enhancing_border
# 3 Macrophage__Malignant     EREG__EGFR:ERBB2 Cell2Cell 6514_enhancing_border
# 4 Macrophage__Malignant    EREG__ERBB2:ERBB4 Cell2Cell 6514_enhancing_border
# 5 Macrophage__Malignant NFASC__CNTN1:CNTNAP1 Cell2Cell 6514_enhancing_border
# 6 Macrophage__Malignant    THY1__ITGAM:ITGB2 Cell2Cell 6514_enhancing_border
obj_cpdb <- obj_cpdb %>%
    select(all_of(common_cols)) %>%
    filter(!is.na(pval)) %>%
    filter(pval < args$alpha) %>%
    select(-pval)
# r$> head(obj_cpdb)
#            source_target complex_interaction        method                Sample
# 1 Macrophage__Macrophage        TGFB1__ITGB1 CellPhoneDBv5 6514_enhancing_border
# 2 Macrophage__Macrophage         SPP1__ITGB1 CellPhoneDBv5 6514_enhancing_border
# 3 Macrophage__Macrophage         CD14__ITGB1 CellPhoneDBv5 6514_enhancing_border
# 4 Macrophage__Macrophage       ADAM17__ITGB1 CellPhoneDBv5 6514_enhancing_border
# 5 Macrophage__Macrophage        ADAM9__ITGB1 CellPhoneDBv5 6514_enhancing_border
# 6 Macrophage__Macrophage       LGALS1__ITGB1 CellPhoneDBv5 6514_enhancing_border

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
# r$> head(interactions_signif)
#                 source_target complex_interaction     method                Sample
# 1        Malignant__Malignant COL4A2__ITGAV:ITGB8 CellChatv2 6514_enhancing_border
# 2  Oligodendrocyte__Malignant COL4A5__ITGAV:ITGB8 CellChatv2 6514_enhancing_border
# 3      Macrophage__Macrophage   SPP1__ITGAV:ITGB1 CellChatv2 6514_enhancing_border
# 4       Microglia__Macrophage   SPP1__ITGAV:ITGB1 CellChatv2 6514_enhancing_border
# 5 Oligodendrocyte__Macrophage   SPP1__ITGAV:ITGB1 CellChatv2 6514_enhancing_border
# 6        Malignant__Malignant   GLS:SLC1A2__GRIA2 CellChatv2 6514_enhancing_border

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
# # A tibble: 6 × 8
# # Groups:   Sample, source_target, complex_interaction [6]
#   Sample                source_target          complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb
#   <chr>                 <chr>                  <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl>
# 1 6514_enhancing_border Macrophage__Macrophage A2M__LRP1                   3        0           1            1       1
# 2 6514_enhancing_border Macrophage__Macrophage ACTR2__ADRB2                1        0           0            1       0
# 3 6514_enhancing_border Macrophage__Macrophage ADAM10__AXL                 1        0           0            1       0
# 4 6514_enhancing_border Macrophage__Macrophage ADAM10__CADM1               1        0           0            1       0
# 5 6514_enhancing_border Macrophage__Macrophage ADAM10__CD44                2        0           0            1       1
# 6 6514_enhancing_border Macrophage__Macrophage ADAM10__GPNMB               3        1           0            1       1

log_info("Take consensus...")
# TODO we could change this if necessary
interactions_mvoted <- all_votes %>%
    group_by(Sample, source_target, complex_interaction) %>%
    mutate(
        lenient_voting = (n_methods >= 3) & (in_liana == 1),
        stringent_voting = (n_methods == 4)
    )
# r$> head(interactions_mvoted)
# # A tibble: 6 × 10
# # Groups:   Sample, source_target, complex_interaction [6]
#   Sample                source_target          complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb lenient_voting stringent_voting
#   <chr>                 <chr>                  <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl> <lgl>          <lgl>
# 1 6514_enhancing_border Macrophage__Macrophage A2M__LRP1                   3        0           1            1       1 FALSE          FALSE
# 2 6514_enhancing_border Macrophage__Macrophage ACTR2__ADRB2                1        0           0            1       0 FALSE          FALSE
# 3 6514_enhancing_border Macrophage__Macrophage ADAM10__AXL                 1        0           0            1       0 FALSE          FALSE
# 4 6514_enhancing_border Macrophage__Macrophage ADAM10__CADM1               1        0           0            1       0 FALSE          FALSE
# 5 6514_enhancing_border Macrophage__Macrophage ADAM10__CD44                2        0           0            1       1 FALSE          FALSE
# 6 6514_enhancing_border Macrophage__Macrophage ADAM10__GPNMB               3        1           0            1       1 TRUE           FALSE

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

saveRDS(
    interactions_df_aggregated,
    glue("{args$output_dir}/{args$sample_id}__interactions_agg_rank.rds")
)

log_info("COMPLETED!")
