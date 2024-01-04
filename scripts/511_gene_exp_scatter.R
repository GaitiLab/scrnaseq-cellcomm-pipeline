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
        description = "Create scatterplots",
    )
    parser$add_argument("--interactions_db",
        type = "character", default = NULL,
        help = "Path to interactions database"
    )
    parser$add_argument("--cutoff_quantile",
        type = "numeric", default = 0.90,
        help = "Cutoff quantile for filtering"
    )
    parser$add_argument("--min_pct_exp",
        type = "numeric", default = 10,
        help = "Minimum percentage of cells expressing a gene"
    )
    parser$add_argument("--res_threshold",
        type = "numeric", default = 1.96,
        help = "Residual threshold for labeling"
    )
    parser$add_argument("--gene_exp_dir",
        type = "character", default = NULL,
        help = "Path to gene expression directory"
    )
    parser$add_argument("--meta",
        type = "character", default = NULL,
        help = "Path to metadata file"
    )
    parser$add_argument("--interactions",
        type = "character", default = NULL,
        help = "Path to interactions file"
    )

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L1"
    args$output_dir <- glue("{here::here()}/output/{args$annot}/511_gene_exp_scatter")
    args$interactions_db <- "001_data_local/interactions_db_v2/ref_db.rds"
    args$cutoff_quantile <- 0.90
    args$min_pct_exp <- 10
    args$res_threshold <- 1.96
    args$gene_exp_dir <- glue("{here::here()}/output/{args$annot}/510_compute_avg_expr")
    args$meta <- glue("{here::here()}/output/{args$annot}/000_data/gbm_regional_study__metadata.rds")
    args$interactions <- glue("{here::here()}/output/{args$annot}/400_consensus/400_samples_interactions_mvoted_w_filters.rds")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)
options(ggrepel.max.overlaps = Inf)

# Load additional libraries
pacman::p_load(Seurat, ggplot2, ggrepel, ggpubr)

colors <- c("sender" = "#0043fc", "receiver" = "#FF3333", "background" = "grey")


log_info("Reading metadata...")
meta <- readRDS(args$meta) %>%
    select(Sample, Region_Grouped) %>%
    distinct()
rownames(meta) <- NULL

log_info("Load detected interactions...")
interactions <- readRDS(args$interactions) %>% separate(source_target, into = c("source", "target"), sep = "__", remove = FALSE) %>% separate(complex_interaction, into = c("ligand_complex", "receptor_complex"), sep = "__", remove = FALSE) %>% rowwise() %>%
mutate(ligand = str_split(ligand_complex, "\\:", simplify = TRUE)[1], receptor = str_split(receptor_complex, "\\:", simplify = TRUE)[1]) %>% ungroup()

# Ligands/receptors from database for shape of points
log_info("Load database of interactions...")
interactions_db <- readRDS(args$interactions_db) %>% rowwise() %>% mutate(ligand = str_split(ligand_complex, "\\:", simplify = TRUE)[1], receptor = str_split(receptor_complex, "\\:", simplify = TRUE)[1]) %>% ungroup()

ligands <- interactions_db %>%
    pull(ligand) %>%
    unique()
receptors <- interactions_db %>%
    pull(receptor) %>%
    unique()

log_info("Load gene expression...")
gene_exp_samples <- list.files(args$gene_exp_dir, full.names = TRUE)
log_info(glue("Found {length(gene_exp_samples)} samples..."))
avg_gene_exp <- do.call(rbind, lapply(gene_exp_samples, function(filename) {
    sample_name <- get_name(filename)
    obj <- readRDS(filename) %>% mutate(Sample = sample_name)
    return(obj)
}))

log_info("Update metadata + filtering on min. pct. cells...")
avg_gene_exp <- avg_gene_exp %>%
    rename(gene = features.plot, cell_type = id) %>%
    filter(pct.exp >= args$min_pct_exp) %>%
    left_join(meta, by = "Sample")

# Average samples per region
overall_avg_gene_exp <- avg_gene_exp %>%
    group_by(Region_Grouped, cell_type, gene) %>%
    summarise(sd = sd(avg.exp), mean = mean(avg.exp)) %>%
    mutate(Region_Grouped = factor(Region_Grouped, levels = c("PT", "TE", "SC")))

# source_oi <- "Malignant"
# target_oi <- "Neuron"
# region_oi <- "PT"

all_regions <- interactions %>%
    filter(stringent_region_pair) %>%
    pull(Region_Grouped) %>%
    unique()

# Count the number of interactions between cell type groups
# Types of filters:
# - stringent_region (voting method based stringent + take into account only region)
# - stringent_region_pair (voting stringent + take into account both region and presence of pair in sample of that region)
# - lenient_region (take into account only region)
# - lenient_region_pair (take into account both region and presence of pair in sample of that region)
options <- c("stringent_region", "stringent_region_pair", "lenient_region", "lenient_region_pair")
# region_oi <- "PT"
for (region_oi in all_regions) {
    available_sender_receiver_pairs <- interactions %>%
        filter(Region_Grouped == region_oi, stringent_region_pair) %>%
        select(source_target, source, target) %>%
        distinct()
    # i <- 1
    for (i in seq_len(nrow(available_sender_receiver_pairs))) {
        source_oi <- available_sender_receiver_pairs$source[i]
        target_oi <- available_sender_receiver_pairs$target[i]
        log_info(glue("Reference: {source_oi} - {target_oi} in {region_oi}..."))
        log_info("Select ligands/receptors from detected interactions...")
        interactions_subset <- interactions %>%
            filter(
                stringent_region_pair,
                Region_Grouped == region_oi, source == source_oi,
                target == target_oi
            ) %>%
            select(ligand, receptor)
        ligands_oi <- interactions_subset %>%
            pull(ligand) %>%
            unique()
        receptors_oi <- interactions_subset %>%
            pull(receptor) %>%
            unique()

        avg_gene_expr_subset <- overall_avg_gene_exp %>%
            filter(cell_type %in% c(source_oi, target_oi)) %>%
            mutate(
                # Ligand or receptor (shape)
                group = ifelse(gene %in% ligands, "ligand", ifelse(gene %in% receptors, "receptor", NA)),
                # Ligand or receptor in list of interest (color)
                exp_by_type =
                    factor(case_when(
                        gene %in% ligands_oi ~ "sender",
                        gene %in% receptors_oi ~ "receiver",
                        TRUE ~ "background"
                    )),
                cell_type = ifelse(cell_type == source_oi, "sender", "receiver")
            )
        if (avg_gene_expr_subset %>% pull(cell_type) %>% unique() %>% length() != 2) {
            log_info("Not enough cell types, skipping...")
            next
        }
        log_info("Convert to wide format...")
        avg_gene_expr_subset_wide <- avg_gene_expr_subset %>%
            select(Region_Grouped, group, exp_by_type, mean, sd, cell_type, gene) %>%
            ungroup() %>%
            pivot_wider(names_from = cell_type, values_from = c(mean, sd), values_fill = NA) %>%
            filter(!is.na(group)) %>%
            mutate(
                sd =
                    case_when(
                        exp_by_type == "sender" ~ sd_sender,
                        exp_by_type == "receiver" ~ sd_receiver,
                        TRUE ~ 1
                    )
            ) %>%
            filter(Region_Grouped == region_oi) %>%
            ungroup() %>%
            select(-Region_Grouped)
        loess_fit <- loess(avg_gene_expr_subset_wide$mean_sender ~ avg_gene_expr_subset_wide$mean_receiver)
        # # Take residuals
        resid <- scale(residuals(loess_fit), scale = TRUE, center = TRUE)
        # his <- hist(resid)

        avg_gene_expr_subset_wide[!(is.na(avg_gene_expr_subset_wide$mean_sender) | is.na(avg_gene_expr_subset_wide$mean_receiver)), "resid"] <- resid
        avg_gene_expr_subset_wide[!(is.na(avg_gene_expr_subset_wide$mean_sender) | is.na(avg_gene_expr_subset_wide$mean_receiver)), "fitted"] <- loess_fit$fitted

        # print(his)
        log_info("Adding labels...")
        avg_gene_expr_subset_wide <- avg_gene_expr_subset_wide %>% mutate(add_label = case_when(
            exp_by_type == "background" ~ NA,
            # abs(sd) < args$res_threshold ~ as.character(gene),
            # abs(sd) < args$res_threshold ~ as.character(gene),
            exp_by_type != "background" ~ as.character(gene),
            TRUE ~ NA
        ))

        log_info("Plotting...")
        plt_all <- ggplot(data = avg_gene_expr_subset_wide) +
            geom_smooth(aes(x = mean_sender, y = mean_receiver),
                se = FALSE, color = "black", linetype = "dashed", method = "loess"
            ) +
            geom_point(aes(
                x = mean_sender, y = mean_receiver, shape = group, color = exp_by_type,
                size = sd,
            )) +
            scale_color_manual(values = colors) +
            custom_theme() +
            labs(
                title = "Expression of ligands & receptors", subtitle = glue("Coloring/labeling based on interactions found between {source_oi} - {target_oi} in {region_oi}"),
                x = glue("log10(Average expression in {source_oi} [sender])"), y = glue("log10(Average expression in {target_oi} [receiver])")
            ) +
            guides(
                color = guide_legend(title = ""), shape = guide_legend(title = ""),
                size = guide_legend(title = "SD")
            ) +
            geom_text_repel(
                aes(mean_sender, mean_receiver,
                    label = add_label, color = exp_by_type
                )
            ) +
            scale_x_log10() +
            scale_y_log10()
        # facet_wrap(~Region_Grouped, ncol = 3)
        plt_all
        log_info("Saving...")
        ggsave(
            plot = plt_all,
            filename = glue("{args$output_dir}/scatter_{source_oi}_{target_oi}_{region_oi}.pdf"),
            width = 8, height = 7
        )
    }
}
