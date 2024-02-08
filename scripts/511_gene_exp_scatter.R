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
    parser$add_argument("--min_pct_exp",
        type = "numeric", default = 10,
        help = "Minimum percentage of cells expressing a gene"
    )
    # parser$add_argument("--res_threshold",
    #     type = "numeric", default = 1.96,
    #     help = "Residual threshold for labeling"
    # )
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
    parser$add_argument("--annot", type = "character", default = "CCI_CellClass_L1", help = "Annotation to use")

    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L1"
    args$output_dir <- glue("{here::here()}/output/{args$annot}/511_gene_exp_scatter")
    args$interactions_db <- glue("{here::here()}/data/interactions_db/ref_db.rds")
    args$min_pct_exp <- 10
    # args$res_threshold <- 1.96
    args$gene_exp_dir <- glue("{here::here()}/output/{args$annot}/510_compute_avg_expr")
    args$meta <- glue("{here::here()}/output/{args$annot}/000_data/gbm_regional_study__metadata.rds")
    args$interactions <- glue("{here::here()}/output/{args$annot}/402_aggregation/402_sample_interactions_mvoted_w_filters.rds")
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
pacman::p_load(Seurat, ggplot2, ggrepel, ggpubr, ggtext)

colors <- c("sender" = "#0043fc", 
"receiver" = "#FF3333", 
"background" = "darkgrey", 
"sender/receiver" = "green")

malignant_id <- ifelse(args$annot == "CCI_CellClass_L1", "Malignant", "_like")

log_info("Reading metadata...")
meta <- readRDS(args$meta) %>%
    select(Sample, Region_Grouped) %>%
    distinct()
rownames(meta) <- NULL
head(meta)
#           Sample Region_Grouped
# 1 6237_2222190_A             PT
# 2 6237_2222190_C             TE
# 3 6237_2222190_D             SC
# 4 6237_2222190_F             SC
# 5 6234_2895153_A             TE
# 6 6234_2895153_B             PT


log_info("Load detected interactions...")
interactions <- readRDS(args$interactions)

# options <- intersect(INTERACTIONS_POST_FILTERING_OPTIONS, colnames(interactions))
# option <- options[1]
interactions <- interactions %>%
    # filter(!!as.symbol(option)) %>%
    separate(source_target, into = c("source", "target"), sep = "__", remove = FALSE) %>%
    separate(complex_interaction, into = c("ligand_complex", "receptor_complex"), sep = "__", remove = FALSE) %>%
    rowwise() %>%
    mutate(ligand = str_split(ligand_complex, "\\:", simplify = TRUE)[1], receptor = str_split(receptor_complex, "\\:", simplify = TRUE)[1]) %>%
    ungroup()
head(interactions)
# # A tibble: 6 × 16
#   complex_interaction ligand_complex receptor_complex Region_Grouped source_target                  source target lenient_voting_samples stringent_voting_sam…¹ lenient_region lenient_region_pair stringent_region stringent_region_pair setname ligand receptor
#   <chr>               <chr>          <chr>            <fct>          <chr>                          <chr>  <chr>  <chr>                  <chr>                  <lgl>          <lgl>               <lgl>            <lgl>                 <chr>   <chr>  <chr>
# 1 A2M__LRP1           A2M            LRP1             PT             Endothelial__Astrocyte         Endot… Astro… 6419_enhancing_border  6419_enhancing_border  FALSE          TRUE                FALSE            TRUE                  default A2M    LRP1
# 2 A2M__LRP1           A2M            LRP1             PT             Endothelial__Differentiated_l… Endot… Diffe… 6419_enhancing_border  6419_enhancing_border  FALSE          TRUE                FALSE            TRUE                  Other-… A2M    LRP1
# 3 A2M__LRP1           A2M            LRP1             PT             Endothelial__Microglia         Endot… Micro… 6419_enhancing_border  6419_enhancing_border  FALSE          TRUE                FALSE            TRUE                  default A2M    LRP1
# 4 A2M__LRP1           A2M            LRP1             PT             Endothelial__Pericyte          Endot… Peric… 6419_enhancing_border  6419_enhancing_border  FALSE          TRUE                FALSE            TRUE                  default A2M    LRP1
# 5 A2M__LRP1           A2M            LRP1             PT             Endothelial__Progenitor_like   Endot… Proge… 6419_enhancing_border  6419_enhancing_border  FALSE          TRUE                FALSE            TRUE                  Other-… A2M    LRP1
# 6 A2M__LRP1           A2M            LRP1             PT             Microglia__Astrocyte           Micro… Astro… 6237_2222190_A, 6419_… 6237_2222190_A, 6419_… FALSE          TRUE                FALSE            TRUE                  default A2M    LRP1
# # ℹ abbreviated name: ¹​stringent_voting_samples

# Ligands/receptors from database for shape of points
log_info("Load database of interactions...")
interactions_db <- readRDS(args$interactions_db) %>%
    rowwise() %>%
    mutate(ligand = str_split(ligand_complex, "\\:", simplify = TRUE)[1], receptor = str_split(receptor_complex, "\\:", simplify = TRUE)[1]) %>%
    ungroup()
head(interactions_db)
# # A tibble: 6 × 9
#   source_genesymbol target_genesymbol interaction  simple_interaction complex_interaction ligand_complex receptor_complex ligand receptor
#   <chr>             <chr>             <chr>        <chr>              <chr>               <chr>          <chr>            <chr>  <chr>
# 1 TNF               TNFRSF1A          TNF_TNFRSF1A TNF__TNFRSF1A      TNF__TNFRSF1A       TNF            TNFRSF1A         TNF    TNFRSF1A
# 2 VIP               VIPR1             VIP_VIPR1    VIP__VIPR1         VIP__VIPR1          VIP            VIPR1            VIP    VIPR1
# 3 HCRT              HCRTR1            HCRT_HCRTR1  HCRT__HCRTR1       HCRT__HCRTR1        HCRT           HCRTR1           HCRT   HCRTR1
# 4 AVP               AVPR1A            AVP_AVPR1A   AVP__AVPR1A        AVP__AVPR1A         AVP            AVPR1A           AVP    AVPR1A
# 5 AVP               AVPR2             AVP_AVPR2    AVP__AVPR2         AVP__AVPR2          AVP            AVPR2            AVP    AVPR2
# 6 CXCL8             CXCR2             CXCL8_CXCR2  CXCL8__CXCR2       CXCL8__CXCR2        CXCL8          CXCR2            CXCL8  CXCR2

ligands <- interactions_db %>%
    pull(ligand) %>%
    unique()
receptors <- interactions_db %>%
    pull(receptor) %>%
    unique()
dual_genes_ref_db <- intersect(ligands, receptors)
log_info(glue("Number of dual genes in reference database: {length(dual_genes_ref_db)}"))

log_info("Load gene expression...")
gene_exp_samples <- list.files(args$gene_exp_dir, full.names = TRUE)
log_info(glue("Found {length(gene_exp_samples)} samples..."))
avg_gene_exp <- do.call(rbind, lapply(gene_exp_samples, function(filename) {
    sample_name <- get_name(filename)
    obj <- readRDS(filename) %>% mutate(Sample = sample_name)
    return(obj)
}))
head(avg_gene_exp)
#      avg.exp    pct.exp features.plot              id avg.exp.scaled         Sample
# 1  9.8958084 60.6986900          NRG2 Progenitor_like    -0.30381134 6234_2895153_A
# 2  0.0000000  0.0000000          EREG Progenitor_like            NaN 6234_2895153_A
# 3  1.2584128 16.5938865        SEMA7A Progenitor_like    -0.70448656 6234_2895153_A
# 4 37.3861813 95.1965066         NFASC Progenitor_like     0.85377580 6234_2895153_A
# 5  1.0427285 11.7903930          TGFA Progenitor_like    -0.42494779 6234_2895153_A
# 6  0.1827182  0.8733624         ANXA1 Progenitor_like    -0.07388671 6234_2895153_A

log_info("Update metadata + filtering on min. pct. cells...")
avg_gene_exp <- avg_gene_exp %>%
    rename(gene = features.plot, cell_type = id) %>%
    # Only take int account expression if gene is expressed in at least X%.
    filter(pct.exp >= args$min_pct_exp) %>%
    left_join(meta, by = "Sample")
head(avg_gene_exp)
#     avg.exp  pct.exp   gene       cell_type avg.exp.scaled         Sample Region_Grouped
# 1  9.895808 60.69869   NRG2 Progenitor_like     -0.3038113 6234_2895153_A             TE
# 2  1.258413 16.59389 SEMA7A Progenitor_like     -0.7044866 6234_2895153_A             TE
# 3 37.386181 95.19651  NFASC Progenitor_like      0.8537758 6234_2895153_A             TE
# 4  1.042729 11.79039   TGFA Progenitor_like     -0.4249478 6234_2895153_A             TE
# 5  5.315499 54.14847   THY1 Progenitor_like      1.1344706 6234_2895153_A             TE
# 6  3.443178 36.68122  L1CAM Progenitor_like     -0.2384047 6234_2895153_A             TE

# Average samples per region
log_info("Compute average expression across samples")
overall_avg_gene_exp <- avg_gene_exp %>%
    group_by(Region_Grouped, cell_type, gene) %>%
    summarise(sd = sd(avg.exp), mean = mean(avg.exp)) %>%
    mutate(Region_Grouped = factor(Region_Grouped, levels = c("PT", "TE", "SC")))
head(overall_avg_gene_exp)
# # A tibble: 6 × 5
# # Groups:   Region_Grouped, cell_type [1]
#   Region_Grouped cell_type gene       sd  mean
#   <fct>          <fct>     <fct>   <dbl> <dbl>
# 1 PT             Microglia NRG2   NA      6.12
# 2 PT             Microglia SEMA7A NA      7.10
# 3 PT             Microglia NFASC   0.851  7.88
# 4 PT             Microglia TGFA    1.90   4.48
# 5 PT             Microglia ANXA1  NA      6.99
# 6 PT             Microglia ICAM1   8.80  12.6

for (option in INTERACTIONS_POST_FILTERING_OPTIONS) {
    if (!option %in% colnames(interactions)) { 
        next
    }
    output_dir <- glue("{args$output_dir}/{option}")
    create_dir(output_dir)
    for (region_oi in REGION_GROUPED_LEVELS) {
        # TODO uncomment for testing
        # region_oi <- REGION_GROUPED_LEVELS[1]
        # i <- 1

        # Only plot the interactions that are detected along the malignant-TME axis (or the other way around)
        available_sender_receiver_pairs <- interactions %>% rowwise() %>% 
            filter(Region_Grouped == region_oi, !!sym(option), str_detect(source_target, malignant_id)) %>%
            select(source_target, source, target) %>%
            distinct()
        #         # A tibble: 6 × 3
        #   source_target              source    target
        #   <chr>                      <chr>     <chr>
        # 1 Astrocyte__Malignant       Astrocyte Malignant
        # 2 Malignant__Astrocyte       Malignant Astrocyte
        # 3 Malignant__Malignant       Malignant Malignant
        # 4 Malignant__Microglia       Malignant Microglia
        # 5 Malignant__Neuron          Malignant Neuron
        # 6 Malignant__Oligodendrocyte Malignant Oligodendrocyte
        if (length(available_sender_receiver_pairs) < 1) {
            next
        }

        for (i in seq_len(nrow(available_sender_receiver_pairs))) {
            # region_oi <- "PT"
            # source_oi <- malignant_id
            # target_oi <- "Neuron"
            source_oi <- available_sender_receiver_pairs$source[i]
            target_oi <- available_sender_receiver_pairs$target[i]
            log_info(glue("Reference: {source_oi} - {target_oi} in {region_oi}..."))
            log_info("Select ligands/receptors from detected interactions...")
            interactions_subset <- interactions %>%
                filter(
                    !!sym(option),
                    Region_Grouped == region_oi, source == source_oi,
                    target == target_oi
                ) %>%
                select(ligand, receptor, source, target, Region_Grouped)
            # # A tibble: 6 × 5
            #   ligand receptor source    target    Region_Grouped
            #   <chr>  <chr>    <chr>     <chr>     <chr>
            # 1 CNTN1  NRCAM    Astrocyte Malignant PT
            # 2 CNTN1  PTPRZ1   Astrocyte Malignant PT
            # 3 FLRT2  ADGRL3   Astrocyte Malignant PT
            # 4 LAMA1  ITGB8    Astrocyte Malignant PT
            # 5 LRIG1  EGFR     Astrocyte Malignant PT
            # 6 NCAM1  PTPRZ1   Astrocyte Malignant PT

            ligands_oi <- interactions_subset %>%
                pull(ligand) %>%
                unique()
            receptors_oi <- interactions_subset %>%
                pull(receptor) %>%
                unique()

            log_info(glue("Number of ligands: {length(ligands_oi)}"))
            log_info(glue("Number of receptors: {length(receptors_oi)}"))

            n_dual_genes <- intersect(ligands_oi, receptors_oi) %>% length()
            log_info(glue("Genes labeled occurring as ligand and receptor: {n_dual_genes}"))

            avg_gene_expr_subset <- overall_avg_gene_exp %>%
                filter(cell_type %in% c(source_oi, target_oi), Region_Grouped == region_oi) %>%
                rowwise() %>%
                mutate(
                    # Ligand or receptor (shape)
                    # group = ifelse(gene %in% ligands, "ligand", ifelse(gene %in% receptors, "receptor", NA)),
                    group = case_when(
                        gene %in% ligands && gene %in% receptors ~ "ligand/receptor",
                        gene %in% ligands ~ "ligand",
                        gene %in% receptors ~ "receptor",
                        TRUE ~ NA
                    ),
                    # Ligand or receptor in list of interest (color)
                    exp_by_type =
                        factor(case_when(
                            gene %in% ligands_oi ~ "sender",
                            gene %in% receptors_oi ~ "receiver",
                            TRUE ~ "background"
                        )),
                    cell_type = ifelse(cell_type == source_oi, "sender", "receiver")
                ) %>%
                ungroup()
            # # A tibble: 6 × 7
            #   Region_Grouped cell_type gene      sd  mean group           exp_by_type
            #   <fct>          <chr>     <fct>  <dbl> <dbl> <chr>           <fct>
            # 1 PT             receiver  NRG2    9.25 13.4  ligand          background
            # 2 PT             receiver  SEMA7A NA     1.94 ligand          background
            # 3 PT             receiver  NFASC  17.1  21.8  ligand          background
            # 4 PT             receiver  ANXA1  NA    18.7  ligand          background
            # 5 PT             receiver  THY1    2.09  4.08 ligand/receptor background
            # 6 PT             receiver  HLA-E  NA     8.54 NA              background

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
                            # If a SD couldn't be computed, use size 1 (otherwise point won't be drawn)
                            TRUE ~ 1
                        )
                ) %>%
                filter(Region_Grouped == region_oi) %>%
                ungroup() %>%
                select(-Region_Grouped)
            # # A tibble: 6 × 8
            #   group           exp_by_type gene   mean_receiver mean_sender sd_receiver sd_sender    sd
            #   <chr>           <fct>       <fct>          <dbl>       <dbl>       <dbl>     <dbl> <dbl>
            # 1 ligand          background  NRG2           13.4         9.84        9.25      4.98     1
            # 2 ligand          background  SEMA7A          1.94       NA          NA        NA        1
            # 3 ligand          background  NFASC          21.8         7.09       17.1       1.98     1
            # 4 ligand          background  ANXA1          18.7         3.19       NA        NA        1
            # 5 ligand/receptor background  THY1            4.08       NA           2.09     NA        1
            # 6 ligand/receptor background  L1CAM           1.59       NA          NA        NA        1

            # Modelling LOESS based on sender/receiver expression
            loess_fit <- loess(avg_gene_expr_subset_wide$mean_sender ~ avg_gene_expr_subset_wide$mean_receiver)

            # # Take residuals
            resid <- scale(residuals(loess_fit), scale = TRUE, center = TRUE)

            avg_gene_expr_subset_wide[!(is.na(avg_gene_expr_subset_wide$mean_sender) | is.na(avg_gene_expr_subset_wide$mean_receiver)), "resid"] <- resid
            avg_gene_expr_subset_wide[!(is.na(avg_gene_expr_subset_wide$mean_sender) | is.na(avg_gene_expr_subset_wide$mean_receiver)), "fitted"] <- loess_fit$fitted
            # # A tibble: 6 × 10
            #   group           exp_by_type gene   mean_receiver mean_sender sd_receiver sd_sender    sd   resid fitted
            #   <chr>           <fct>       <fct>          <dbl>       <dbl>       <dbl>     <dbl> <dbl>   <dbl>  <dbl>
            # 1 ligand          background  NRG2           13.4         9.84        9.25      4.98     1 -0.0907   14.1
            # 2 ligand          background  SEMA7A          1.94       NA          NA        NA        1 NA        NA
            # 3 ligand          background  NFASC          21.8         7.09       17.1       1.98     1 -0.353    22.0
            # 4 ligand          background  ANXA1          18.7         3.19       NA        NA        1 -0.404    20.1
            # 5 ligand/receptor background  THY1            4.08       NA           2.09     NA        1 NA        NA
            # 6 ligand/receptor background  L1CAM           1.59       NA          NA        NA        1 NA        NA

            log_info("Adding labels...")
            avg_gene_expr_subset_wide <- avg_gene_expr_subset_wide %>% mutate(add_label = case_when(
                exp_by_type == "background" ~ NA,
                # abs(sd) < args$res_threshold ~ as.character(gene),
                # abs(sd) < args$res_threshold ~ as.character(gene),
                exp_by_type != "background" ~ as.character(gene),
                TRUE ~ NA
            ))
            # # A tibble: 6 × 11
            #   group    exp_by_type gene  mean_receiver mean_sender sd_receiver sd_sender    sd   resid fitted add_label
            #   <chr>    <fct>       <fct>         <dbl>       <dbl>       <dbl>     <dbl> <dbl>   <dbl>  <dbl> <chr>
            # 1 ligand   background  NRG2          13.4         9.84        9.25      4.98     1 -0.0907   14.1 NA
            # 2 ligand   background  SEMA…          1.94       NA          NA        NA        1 NA        NA   NA
            # 3 ligand   background  NFASC         21.8         7.09       17.1       1.98     1 -0.353    22.0 NA
            # 4 ligand   background  ANXA1         18.7         3.19       NA        NA        1 -0.404    20.1 NA
            # 5 ligand/… background  THY1           4.08       NA           2.09     NA        1 NA        NA   NA
            # 6 ligand/… background  L1CAM          1.59       NA          NA        NA        1 NA        NA   NA

            log_info("Plotting...")
            plt_all <- ggplot(data = avg_gene_expr_subset_wide) +
                geom_smooth(aes(x = mean_sender, y = mean_receiver),
                    se = TRUE, color = "black", linetype = "dashed", method = "loess"
                ) +
                geom_point(aes(
                    # Shape = ligand/receptor based on database
                    # Color = expressed by sender or receiver cell type
                    x = mean_sender, y = mean_receiver, shape = group, color = exp_by_type,
                    size = sd,
                )) +
                scale_color_manual(values = colors) +
                custom_theme() +
                labs(
                    title = "Expression of ligands & receptors",
                    subtitle = glue("Coloring/labeling based on interactions found between {source_oi} - {target_oi} in {region_oi}"),
                    x = glue("Average expression of {source_oi} (log<sub>10</sub>)"),
                    y = glue("Average expression of {target_oi} (log<sub>10</sub>)")
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
                scale_x_log10(
                    breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))
                ) +
                scale_y_log10(
                    breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))
                )

            plt_all
            log_info("Saving...")
            ggsave(
                plot = plt_all,
                filename = glue("{output_dir}/scatter_{source_oi}_{target_oi}_{region_oi}.pdf"),
                width = 10, height = 10
            )
            auto_crop(glue("{output_dir}/scatter_{source_oi}_{target_oi}_{region_oi}.pdf"))
        }
    }
}
