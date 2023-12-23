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
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/CellClass_L4_min3_types/521_gene_exp_scatter")
    args$interactions_db <- "001_data_local/interactions_db/interactions_ref.rds"
    args$stringency <- 0
    args$cutoff_quantile <- 0.90
    args$min_pct_exp <- 10
    args$res_threshold <- 1.96
    args$gene_exp_dir <- glue("{here::here()}/output/CellClass_L4_min3_types/520_compute_avg_expr")
    args$meta <- glue("{here::here()}/001_data_local/seurat_annot_adapted__metadata.rds")
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
    select(Sample, Region) %>%
    distinct()
rownames(meta) <- NULL

log_info("Load detected interactions...")
interactions <- readRDS(glue("{here::here()}/output/CellClass_L4_min3_types/400_consensus/400c_post_filtering.rds"))

# Ligands/receptors from database for shape of points
log_info("Load database of interactions...")
interactions_db <- readRDS(args$interactions_db)
ligands <- interactions_db %>%
    pull(genename_a) %>%
    unique()
receptors <- interactions_db %>%
    pull(genename_b) %>%
    unique()

log_info("Load gene expression...")
gene_exp_samples <- list.files(args$gene_exp_dir, full.names = TRUE)
avg_gene_exp <- do.call(rbind, lapply(gene_exp_samples, function(filename) {
    sample_name <- get_name(filename)
    obj <- readRDS(filename) %>% mutate(Sample = sample_name)
    return(obj)
}))

# Update metadata + filtering on min. pct. cells
avg_gene_exp <- avg_gene_exp %>%
    rename(gene = features.plot, cell_type = id) %>%
    filter(pct.exp >= args$min_pct_exp) %>%
    left_join(meta, by = "Sample")

# Average samples per region
overall_avg_gene_exp <- avg_gene_exp %>%
    group_by(Region, cell_type, gene) %>%
    summarise(sd = sd(avg.exp), mean = mean(avg.exp)) %>%
    mutate(Region = factor(Region, levels = c("PT", "TE", "SC")))

# source_oi <- "Malignant"
# target_oi <- "Neuron"
# region_oi <- "PT"

all_regions <- interactions %>%
    filter(is_stringent == args$stringency) %>%
    pull(Region) %>%
    unique()


for (region_oi in all_regions) {
    available_sender_receiver_pairs <- interactions %>%
        filter(Region == region_oi, is_stringent == args$stringency, cond_min_samples_region) %>%
        select(source_target, source, target) %>%
        distinct()

    for (i in seq_len(nrow(available_sender_receiver_pairs))) {
        source_oi <- available_sender_receiver_pairs$source[i]
        target_oi <- available_sender_receiver_pairs$target[i]
        log_info(glue("Reference: {source_oi} - {target_oi} in {region_oi}..."))
        log_info("Select ligands/receptors from detected interactions...")
        interactions_subset <- interactions %>%
            filter(
                cond_min_samples_region,
                is_stringent == args$stringency,
                Region == region_oi, source == source_oi,
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
            select(Region, group, exp_by_type, mean, sd, cell_type, gene) %>%
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
            filter(Region == region_oi) %>%
            ungroup() %>%
            select(-Region)


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
        # facet_wrap(~Region, ncol = 3)
        plt_all
        log_info("Saving...")
        ggsave(
            plot = plt_all,
            filename = glue("{args$output_dir}/scatter_{source_oi}_{target_oi}_{region_oi}_{args$stringency}.pdf"),
            width = 8, height = 7
        )
    }
}

# plt_all
# plt_sub <- plt_all + xlim(0, 20) + ylim(0, 20) +  geom_text_repel(
#         aes(mean_Malignant, mean_Neuron,
#             # label = ifelse((abs(resid) > args$res_threshold & !is.na(exp_by_type)), as.character(gene), ""),
#             label = ifelse(exp_by_type == "background", "", as.character(gene)), color = exp_by_type
#         ),) + labs(subtitle = "Zoomed in")
# plt_combi <- ggarrange(plt_all, plt_sub, ncol = 2, common.legend = TRUE, legend = "bottom")
# plt_combi
# ggsave(plot = plt_all, filename = glue("{args$output_dir}/scatter_{region_oi}_{args$stringency}.pdf"), width = 10, height = 5)

# # ---- Filter by cutoff quantile ---- #
# max_source <- quantile(avg_expr_subset[[source_oi]], args$cutoff_quantile)
# max_target <- quantile(avg_expr_subset[[target_oi]], args$cutoff_quantile)
# expr_cutoff <- ceiling(max(c(max_source, max_target)))

# avg_expr_subset_filtered <- avg_expr_subset %>% filter(Malignant <= expr_cutoff, Neuron <= expr_cutoff)
# loess_fit <- loess(avg_expr_subset_filtered$Malignant ~ avg_expr_subset_filtered$Neuron)

# # Take residuals
# resid <- scale(residuals(loess_fit), scale = TRUE, center = TRUE)
# his <- hist(resid)
# print(his)

# avg_expr_subset_filtered$resid <- resid
# avg_expr_subset_filtered$fitted <- loess_fit$fitted

# plt_subset <- ggplot(data = avg_expr_subset_filtered) +
#     geom_smooth(aes(x = Malignant, y = Neuron),
#         se = FALSE, color = "lightgrey", linetype = "dashed", method = "loess"
#     ) +
#     geom_point(aes(x = Malignant, y = Neuron, shape = group, color = exp_by_type)) +
#     custom_theme() +
#     labs(x = glue("Average expression in {source_oi} (sender)"), y = glue("Average expression in {target_oi} (receiver)")) +
#     guides(color = guide_legend(title = ""), shape = guide_legend(title = "")) +
#     # Only label the genes with residuals > 1.96 and
#     # Color labels of genes that are also within the ligand/receptor list of interest in black
#     geom_text_repel(
#         aes(Malignant, Neuron,
#             label = ifelse((abs(resid) > args$res_threshold & !is.na(exp_by_type)), as.character(gene), "")
#         )
#     )
# ggsave(plot = plt_subset, filename = glue("{args$output_dir}/scatter_sub_{region_oi}_{args$stringency}.pdf"), width = 10, height = 5)


# ggarrange(plt_all, plt_subset)
# # interactions_with_exp <- interactions %>%
#     left_join(avg_expr %>% rename(source_avg_exp = avg.exp, source_pct_exp = pct.exp, source_avg_exp_scaled = avg.exp.scaled), by = c("source" = "id", "ligand" = "gene")) %>%
#     left_join(avg_expr %>% rename(target_avg_exp = avg.exp, target_pct_exp = pct.exp, target_avg_exp_scaled = avg.exp.scaled), by = c("target" = "id", "receptor" = "gene"))

# interactions_with_exp_group <- interactions_with_exp %>% filter(is_stringent == 0, cond_min_samples_region, Region == "PT", source == "Malignant", target == "Neuron")


# loess_fit <- loess(interactions_with_exp_group$source_avg_exp ~ interactions_with_exp_group$target_avg_exp)

# # Take residuals
# resid <- scale(residuals(loess_fit), scale = TRUE, center = TRUE)
# # his <- histogram(resid)
# # print(his)

# interactions_with_exp_group$resid <- resid
# interactions_with_exp_group$fitted <- loess_fit$fitted
# interactions_with_exp_group$type <- ""
# interactions_with_exp_group[interactions_with_exp_group$gene %in% interactinos_db$genename_a, "type"] <- "ligands"
# interactions_with_exp_group[interactions_with_exp_group$gene %in% interactions_db$genename_b, "type"] <- "receptors"



# p <- ggplot(data = avg_expr, aes(source_avg_exp, target_avg_exp)) +
#     geom_smooth(aes(source_avg_exp, target_avg_exp),
#         se = FALSE, color = "lightgrey", linetype = "dashed", method = "loess"
#     ) +
#     geom_point(aes(color = resid, shape = type)) +
#     scale_color_gradient2(low = "#0043fc", mid = "white", high = "#FF3333") +
#     guides(color = guide_colorbar("Residual")) +
#     labs(title = glue("Average expression of ligand/receptor genes in \n{str_replace(source_target, '__', '-')} interactions"), subtitle = glue("{str_replace_all(region, '_', ' ')}"), x = glue("{c(str_split(source_target, '__', simplify = TRUE)[1])} cell gene expression"), y = glue("{c(str_split(source_target, '__', simplify = TRUE)[2])} cell gene expression")) +
#     scale_x_continuous(limits = c(0, 3)) +
#     scale_y_continuous(limits = c(0, 3)) +
#     theme(
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")
#     ) +
#     # Only label the genes with residuals > 1.96 and
#     # Color labels of genes that are also within the ligand/receptor list of interest in black
#     geom_text_repel(
#         aes(source_avg_exp, target_avg_exp,
#             label = ifelse((abs(resid) > args$res_threshold), as.character(gene), "")
#         ),
#         color = ifelse(((abs(resid) > args$res_threshold) & (interactions_with_exp_group$gene %in% genes_oi$gene)), "black", "grey")
#     )
