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
        description = "Visualizing average expression of ligands/receptors",
    )
    parser$add_argument("--gene_expr", type = "character", help = "Path to gene expression matrix")
    parser$add_argument("--interactions_ref", type = "character", help = "Path to interactions reference")
    parser$add_argument("--interactions_oi", type = "character", help = "Path to interactions of interest")
    parser$add_argument("--min_z_score", type = "numeric", help = "Minimum z-score to consider a gene as highly expressed")
    parser$add_argument("--min_pct_expressed", type = "numeric", help = "Minimum percentage of cells expressing a gene")

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/final_output/CellClass_L2_all/")
    args$min_z_score <- 0.5
    args$min_pct_expressed <- 0
    args$is_stringent <- 1
    # args$samples_oi <-
    # glue("{here::here()}/final_output/CellClass_L2_all/samples_oi_min200.txt")
    args$sample_oi <- NULL
    args$input_dir <- glue("{here::here()}/final_output/CellClass_L2_all/501a_compute_avg_expr")
    args$interactions_ref <- glue("{here::here()}/data/interactions_db/interactions_ref.rds")
    args$interactions_oi <- glue("{here::here()}/final_output/CellClass_L2_all/402_post_filtering/ccis_post_filtering__stringency_{args$is_stringent}_groupby_1.rds")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- paste0(args$output_dir, "/501b_lineplots")
create_dir(output_dir)

# Additional libraries
pacman::p_load(tidyr, Seurat, dplyr, ggplot2, ggrepel, scales, ggpubr, grid, pbapply)

#  TODO remove once using new objects
rename_dict <- c("AC-like" = "Malignant_AC", "MES-like" = "Malignant_MES", "NPC & OPC-like" = "Malignant_NPC_&_OPC", "TAM-BDM" = "Macrophage")

all_interactions <- readRDS(args$interactions_oi)

regions_oi <- all_interactions %>%
    pull(Region) %>%
    unique()


for (region in regions_oi) {
    interactions <- all_interactions %>% filter(Region == region)

    log_info("Add the main ligand and receptor for each interactions...")
    interactions_db <- readRDS(args$interactions_ref) %>% mutate(interaction = str_replace_all(interaction, "__", " - "))
    main_lr <- merge(interactions, interactions_db, by = c("interaction"), all.x = TRUE)
    main_lr <- main_lr %>% distinct(, .keep_all = TRUE)
    log_info(glue("Number of interactions: {nrow(main_lr)}"))

    source_targets_oi <- interactions %>%
        pull(source_target) %>%
        unique()

    log_info("Load average expression...")
    if (!is.null(args$samples_oi)) {
        samples_oi <- read.table(args$samples_oi) %>% pull(V1)
        filepaths <- paste0(args$input_dir, "/", samples_oi, "__lr_avg_expr.rds")
    } else {
        filepaths <- list.files(args$input_dir, full.names = TRUE, pattern = glue("*__lr_avg_expr.rds"))
    }
    lr_avg_expr <- do.call(rbind, lapply(seq_along(filepaths), function(i) {
        obj <- readRDS(filepaths[i]) %>% mutate(
            sample_id = samples_oi[i],
            # custom_annot = str_replace_all(custom_annot, rename_dict))
        return(obj)
    }))


    gene_expr_avg_long <- lr_avg_expr %>%
        group_by(CellClass_L2, gene) %>%
        summarise(avg_expr = mean(avg_expr))

    log_info("Only keep interactions for which both ligand and receptor are expressed...")
    mask <- sapply(seq_len(nrow(main_lr)), function(i) {
        ligand <- main_lr[i, "ligand"]
        receptor <- main_lr[i, "receptor"]
        source <- main_lr[i, "source"]
        target <- main_lr[i, "target"]

        ligand_expr <- nrow(gene_expr_avg_long %>% filter(custom_annot == source, gene == ligand)) > 0

        receptor_expr <- nrow(gene_expr_avg_long %>% filter(custom_annot == target, gene == receptor)) > 0

        return(ligand_expr & receptor_expr)
    })

    main_lr <- main_lr[mask, ]

    log_info("Merge average expression with main ligand-receptor...")
    main_lr$pair_id <- seq_len(nrow(main_lr))
    lr_avg_expr <- merge(main_lr, gene_expr_avg_long,
        by.x = c("genename_a", "source"), by.y = c("gene", "custom_annot"),
        all.x = TRUE
    )
    head(lr_avg_expr)
    log_debug("Ligand average expression...")
    lr_avg_expr <- lr_avg_expr %>% rename(source_avg_expr_ligand = avg_expr)
    lr_avg_expr <- merge(lr_avg_expr,
        gene_expr_avg_long,
        by.x = c("genename_a", "target"),
        by.y = c("gene", "custom_annot"), all.x = TRUE
    )
    lr_avg_expr <- lr_avg_expr %>% rename(target_avg_expr_ligand = avg_expr)

    log_debug("Receptor average expression...")
    lr_avg_expr <- merge(lr_avg_expr, gene_expr_avg_long,
        by.x = c("genename_b", "source"),
        by.y = c("gene", "custom_annot"), all.x = TRUE
    )
    lr_avg_expr <- lr_avg_expr %>% rename(source_avg_expr_receptor = avg_expr)
    lr_avg_expr <- merge(lr_avg_expr, gene_expr_avg_long,
        by.x = c("genename_b", "target"), by.y = c("gene", "custom_annot"), all.x = TRUE
    )
    lr_avg_expr <- lr_avg_expr %>%
        rename(target_avg_expr_receptor = avg_expr) %>%
        mutate(
            source_avg_expr_ligand = case_when(
                abs(source_avg_expr_ligand) > 3 & !is.na(source_avg_expr_ligand) ~ sign(source_avg_expr_ligand * 3), TRUE ~ source_avg_expr_ligand
            ),
            target_avg_expr_ligand = case_when(
                abs(target_avg_expr_ligand) > 3 & !is.na(target_avg_expr_ligand) ~ sign(target_avg_expr_ligand * 3), TRUE ~ target_avg_expr_ligand
            ),
            source_avg_expr_receptor = case_when(
                abs(source_avg_expr_receptor) > 3 & is.na(source_avg_expr_receptor) ~ sign(source_avg_expr_receptor * 3), TRUE ~ source_avg_expr_receptor
            ),
            target_avg_expr_receptor = case_when(
                abs(target_avg_expr_receptor) > 3 & !is.na(target_avg_expr_receptor) ~ sign(target_avg_expr_receptor * 3), TRUE ~ target_avg_expr_receptor
            ),
            pair_id = seq_len(nrow(.))
        )

    log_info("Visualize expression patterns...")
    list_of_plots <- list()
    all_combis <- list()

    for (source_target in source_targets_oi) {
        log_info(glue("Select interactions for {source_target}..."))
        source <- str_split(source_target, "__", simplify = TRUE)[1]
        target <- str_split(source_target, "__", simplify = TRUE)[2]
        source_label <- glue("{source}\n ligand")
        target_label <- glue("{target}\n receptor")
        tmp_sub <- lr_avg_expr[lr_avg_expr$source_target == source_target, ]

        log_info("Format data...")
        log_info("Format receptor data...")
        tmp_receptor <- tmp_sub %>%
            select(pair_id, target_avg_expr_receptor, genename_b, source_target, setname) %>%
            rename(avg_expr = target_avg_expr_receptor, genename = genename_b)
        tmp_receptor$type <- target_label
        tmp_receptor$genename <- as.factor(tmp_receptor$genename)
        tmp_receptor$gene_id <- as.numeric(unclass(tmp_receptor$genename))
        tmp_receptor$gene_id <- tmp_receptor$gene_id

        log_info("Format ligand data...")
        tmp_ligand <- tmp_sub %>%
            select(pair_id, source_avg_expr_ligand, genename_a, source_target, setname) %>%
            rename(avg_expr = source_avg_expr_ligand, genename = genename_a)
        tmp_ligand$type <- source_label
        tmp_ligand$genename <- as.factor(tmp_ligand$genename)
        tmp_ligand$gene_id <- as.numeric(unclass(tmp_ligand$genename))
        log_info("Combine data...")
        tmp_combi <- rbind(tmp_ligand, tmp_receptor)

        log_info("Format data...")
        tmp_combi$type <- as.factor(tmp_combi$type)
        tmp_combi$genename <- as.factor(tmp_combi$genename)

        log_debug("Control alignment/position of labels...")
        tmp_combi[tmp_combi$type == source_label, "nudge_x"] <- -0.1
        tmp_combi[tmp_combi$type == target_label, "nudge_x"] <- 0.1
        tmp_combi[tmp_combi$type == source_label, "hjust"] <- "right"
        tmp_combi[tmp_combi$type == target_label, "hjust"] <- "left"
        tmp_combi$type <- factor(tmp_combi$type, levels = c(source_label, target_label))
        all_combis[[source_target]] <- tmp_combi
    }

    all_combis <- do.call("rbind", all_combis)
    rownames(all_combis) <- NULL

    avail_celltypes <- lr_avg_expr %>%
        select(source, target) %>%
        pull() %>%
        unique()
    celltype_levels <- c(paste0(avail_celltypes, "\nligand"), paste0(avail_celltypes, "\nreceptor"))

    all_combis <- all_combis %>% mutate(
        type = factor(type, levels = celltype_levels),
        # For binarizing
        dummy = (avg_expr >= args$min_z_score) * 1,
    )

    # Only keep the pairs of which ligand and receptor are highly expressed (>= threshold)
    ids <- all_combis %>%
        select(pair_id, dummy) %>%
        group_by(pair_id) %>%
        summarise(sum = sum(dummy)) %>%
        filter(sum == 2) %>%
        pull(pair_id)

    all_combis <- all_combis %>% mutate(is_highly_expressed = ifelse(pair_id %in% ids, "yes", "no"))

    all_combis %>%
        group_by(source_target, pair_id) %>%
        reframe(interaction = paste0(genename, collapse = " - "))

    log_info("Visualize")
    p <- create_lineplot(all_combis)

    max_n <- all_combis %>%
        select(source_target, pair_id) %>%
        distinct() %>%
        group_by(source_target) %>%
        summarise(n = n()) %>%
        pull(n) %>%
        max()

    ggsave(glue("{output_dir}/overall_lineplot_lr__{region}__stringency_{args$is_stringent}.pdf"), plot = p, width = 14, height = round(max_n * .35))
}
