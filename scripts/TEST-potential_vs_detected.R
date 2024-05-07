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
        description = "Inspect relationship between potential and detected number of interactions"
    )
    parser$add_argument("--interactions_db", type = "character", help = "INteractions database")
    parser$add_argument("--meta", type = "character", help = "Path to metadata RDS file")
    parser$add_argument("--interactions", type = "character", help = "Path to interactions mvoted RDS file")
    parser$add_argument("--gene_exp_dir", type = "character", help = "Path to directory with average gene expression for each sample (list of RDS files)")
    parser$add_argument("--annot", type = "character", help = "Annotatation used")
    parser$add_argument("--min_pct", type = "numeric", help = "Minimum percentage of cells expressing a gene", default = 0.10)
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L1"
    args$output_dir <- glue("{here::here()}/output/{args$annot}/TEST-potential_vs_detected")
    args$meta <- glue("{here::here()}/output/CCI_CellClass_L2/000_data/bw_gbm_regional_study__metadata.rds")
    args$min_pct <- 0
    args$interactions_db <- glue("{here::here()}/data/interactions_db/ref_db.rds")
    args$interactions <- glue("{here::here()}/output/{args$annot}/401_combine_samples/401_samples_interactions_mvoted.rds")

    args$gene_exp_dir <- glue("{here::here()}/output/{args$annot}/510_compute_avg_expr")
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
pacman::p_load(Seurat, ggplot2, ggpubr, ggtext)

log_info("Load interactions database...")
interactions_db <- readRDS(args$interactions_db) %>% separate(simple_interaction, into = c("ligand", "receptor"), sep = "__", remove = FALSE)

log_info("Load interactions...")
interactions <- readRDS(args$interactions) %>% separate(source_target, into = c("source", "target"), sep = "__", remove = FALSE)
# # A tibble: 6 × 16
#   Sample         source_target        source    target    complex_interaction n_methods in_liana in_cellchat in_cell2cell in_cpdb lenient_voting stringent_voting Patient Region_Grouped Batch   Platform
#   <chr>          <chr>                <chr>     <chr>     <chr>                   <int>    <dbl>       <dbl>        <dbl>   <dbl> <lgl>          <lgl>            <chr>   <chr>          <chr>   <chr>
# 1 6234_2895153_A Malignant__Malignant Malignant Malignant A2M__LRP1                   1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio
# 2 6234_2895153_A Malignant__Malignant Malignant Malignant ACE__BDKRB2                 1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio
# 3 6234_2895153_A Malignant__Malignant Malignant Malignant ADAM10__AXL                 1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio
# 4 6234_2895153_A Malignant__Malignant Malignant Malignant ADAM10__CADM1               1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio
# 5 6234_2895153_A Malignant__Malignant Malignant Malignant ADAM10__GPNMB               1        0           0            1       0 FALSE          FALSE            6234    TE             Batch_2 ParseBio

log_info("Load metadata...")
meta <- readRDS(args$meta)
rownames(meta) <- NULL
meta <- meta %>%
    select(Sample, Region_Grouped) %>%
    distinct()

log_info("Load database with interactions...")
db <- readRDS(args$interactions_db)

genes_db <- db %>%
    select(source_genesymbol, target_genesymbol) %>%
    pull() %>%
    unique() %>%
    str_split("_", simplify = TRUE) %>%
    unlist()
genes_db <- unique(genes_db[genes_db != ""])
log_info(glue("Number of genes: {length(genes_db)}"))

log_info("Compute number of interactions per pair per sample...")
types_of_voting <- c("lenient_voting", "stringent_voting")

for (type_of_voting in types_of_voting) {
    # type_of_voting <- types_of_voting[1]
    all_n_interactions_by_sample <- interactions %>%
        select(Sample, source_target, source, target, complex_interaction, Region_Grouped, !!sym(type_of_voting)) %>%
        distinct() %>%
        filter(!!as.symbol(type_of_voting)) %>%
        group_by(Region_Grouped, Sample, source_target, source, target) %>%
        reframe(n_interactions = n())
    # # A tibble: 6 × 6
    #   Region_Grouped Sample         source_target              source    target          n_interactions
    #   <chr>          <chr>          <chr>                      <chr>     <chr>                    <int>
    # 1 PT             6234_2895153_B Malignant__Malignant       Malignant Malignant                  156
    # 2 PT             6234_2895153_B Malignant__Microglia       Malignant Microglia                  174
    # 3 PT             6234_2895153_B Malignant__Neuron          Malignant Neuron                      89
    # 4 PT             6234_2895153_B Malignant__Oligodendrocyte Malignant Oligodendrocyte             97
    # 5 PT             6234_2895153_B Microglia__Malignant       Microglia Malignant                  143
    # 6 PT             6234_2895153_B Microglia__Microglia       Microglia Microglia                  258

    files_avg_gene_expr_sample <- list.files(args$gene_exp_dir, full.names = TRUE)
    avg_gene_expr <- lapply(files_avg_gene_expr_sample, function(file) {
        avg_gene_exp <- readRDS(file)
        number_of_genes_per_class <- avg_gene_exp %>%
            filter(pct.exp > args$min_pct, features.plot %in% genes_db) %>%
            mutate(dummy = 1) %>%
            group_by(id) %>%
            reframe(n = sum(dummy), expressed_genes = paste0(features.plot, collapse = ", ")) %>%
            mutate(fraction = n / length(genes_db), Sample = get_name(file))
    })
    avg_gene_expr <- do.call(rbind, avg_gene_expr)
    # # A tibble: 6 × 5
    #   id            n expressed_genes                                                                                                                                                                                                 fraction Sample
    #   <fct>     <dbl> <chr>                                                                                                                                                                                                              <dbl> <chr>
    # 1 Malignant   842 THY1, GP1BA, CDH1, L1CAM, JAM2, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL18, NT5E, AGRN, APP, BMPR1B, BMPR2, FLRT1, FLRT3, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITGB3, ITGB5, B2M, NRXN1, CD47, CNTN2, ITGA4, ITGA…    0.714 6234_…
    # 2 Microglia   793 THY1, L1CAM, JAM2, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL18, NT5E, AGRN, APP, BMPR1B, BMPR2, FLRT1, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITGB3, ITGB5, B2M, NRXN1, CD47, CNTN2, ITGA4, ITGA9, PODXL, ADGRB1, AD…    0.673 6234_…
    # 3 Neuron      956 THY1, GP1BA, CDH1, VCAM1, L1CAM, JAM2, MADCAM1, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL18, CRLF2, NT5E, SLC18A2, AGRN, APP, BMPR1B, BMPR2, FLRT1, FLRT3, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITGB3, ITGB5, B2M,…    0.811 6234_…
    # 4 Malignant   912 THY1, GP1BA, CDH1, VCAM1, L1CAM, JAM2, MADCAM1, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL18, TSLP, NT5E, SLC18A2, AGRN, APP, BMPR1B, BMPR2, FLRT1, FLRT3, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITGB3, ITGB5, B2M, …    0.774 6234_…
    # 5 Microglia   831 THY1, CDH1, L1CAM, JAM2, MADCAM1, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL18, NT5E, SLC18A2, AGRN, APP, BMPR1B, BMPR2, FLRT1, FLRT3, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITGB3, ITGB5, B2M, NRXN1, CD47, CNTN2, …    0.705 6234_…
    # 6 Neuron     1001 THY1, FCER2, GP1BA, CDH1, VCAM1, L1CAM, JAM2, MADCAM1, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL13, IL18, NT5E, SLC18A2, AGRN, APP, BMPR1B, BMPR2, FLRT1, FLRT3, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITGB3, ITGB5…    0.849 6234_…

    all_df <- avg_gene_expr %>%
        left_join(meta %>% distinct(), by = "Sample")
    # # A tibble: 6 × 6
    #   id            n expressed_genes                                                                                                                                                                                  fraction Sample Region_Grouped
    #   <fct>     <dbl> <chr>                                                                                                                                                                                               <dbl> <chr>  <chr>
    # 1 Malignant   842 THY1, GP1BA, CDH1, L1CAM, JAM2, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL18, NT5E, AGRN, APP, BMPR1B, BMPR2, FLRT1, FLRT3, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITGB3, ITGB5, B2M, NRXN1, CD47, CNT…    0.714 6234_… TE
    # 2 Microglia   793 THY1, L1CAM, JAM2, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL18, NT5E, AGRN, APP, BMPR1B, BMPR2, FLRT1, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITGB3, ITGB5, B2M, NRXN1, CD47, CNTN2, ITGA4, ITGA9, PO…    0.673 6234_… TE
    # 3 Neuron      956 THY1, GP1BA, CDH1, VCAM1, L1CAM, JAM2, MADCAM1, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL18, CRLF2, NT5E, SLC18A2, AGRN, APP, BMPR1B, BMPR2, FLRT1, FLRT3, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITG…    0.811 6234_… TE
    # 4 Malignant   912 THY1, GP1BA, CDH1, VCAM1, L1CAM, JAM2, MADCAM1, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL18, TSLP, NT5E, SLC18A2, AGRN, APP, BMPR1B, BMPR2, FLRT1, FLRT3, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITGB…    0.774 6234_… PT
    # 5 Microglia   831 THY1, CDH1, L1CAM, JAM2, MADCAM1, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL18, NT5E, SLC18A2, AGRN, APP, BMPR1B, BMPR2, FLRT1, FLRT3, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB1, ITGB3, ITGB5, B2M, NRXN1…    0.705 6234_… PT
    # 6 Neuron     1001 THY1, FCER2, GP1BA, CDH1, VCAM1, L1CAM, JAM2, MADCAM1, COL13A1, CSPG4, JAM3, F11R, CRLF1, IL13, IL18, NT5E, SLC18A2, AGRN, APP, BMPR1B, BMPR2, FLRT1, FLRT3, SLC17A7, ALOX5, PLAUR, ITGAV, ITGB…    0.849 6234_… PT
    saveRDS(all_df, glue("{args$output_dir}/genes_per_class_{type_of_voting}.rds"))


    # log_info(glue("Number of samples: {length(all_df %>% pull(Sample) %>% unique())}"))

    # p_fraction <- ggplot(all_df, aes(x = id, y = fraction, color = Region_Grouped)) +
    #     geom_boxplot(outlier.shape = NA) +
    #     coord_flip() +
    #     default_theme() +
    #     labs(x = "Cell type", title = "Number of ligands/receptors expressed from database per cell type", subtitle = glue("min. fraction of cells that express the ligand/receptor > {args$min_pct}"), y = "Fraction of ligands and receptors expressed") +
    #     geom_point(position = position_jitterdodge(), size = .75)

    # p_fraction_stat_test <- ggboxplot(all_df, x = "id", y = "fraction",
    # facet.by = "Region_Grouped", add = "jitter", outlier.shape = NA, ylab = "Fraction of ligands and receptors expressed",
    # xlab = "Cell type", title = "Number of ligands/receptors expressed from database per cell type",
    # subtitle = glue("min. fraction of cells that express the ligand/receptor > {args$min_pct}"), legend = "bottom", ggtheme = theme_pubr(), ncol = 5) +
    # stat_compare_means() + default_theme() + theme(axis.text.x)
    # p_fraction_stat_test

    # ggsave(plot = p_fraction_stat_test, filename = glue("{args$output_dir}/fractiongenes_per_class_minpct{args$min_pct * 100}__stat_test.pdf"), width = 10, height = 10, dpi = 300)


    # p_fraction
    # p_counts <- ggplot(all_df, aes(x = id, y = n, color = Region_Grouped)) +
    #     geom_boxplot(outlier.shape = NA) +
    #     coord_flip() +
    #     default_theme() +
    #     labs(x = "Cell type", title = "Number of ligands/receptors expressed from database per cell type", subtitle = glue("min. fraction of cells that express the ligand/receptor > {args$min_pct}"), y = "Number of genes") +
    #     geom_point(position = position_jitterdodge(), size = .75)

    # ggsave(plot = p_fraction, filename = glue("{args$output_dir}/fractiongenes_per_class_minpct{args$min_pct * 100}.pdf"), width = 10, height = 12, dpi = 300)

    # ggsave(plot = p_counts, filename =
    # glue("{args$output_dir}/ngenes_per_class_minpct{args$min_pct * 100}.pdf"),
    # width = 10, height = 12, dpi = 300)

    # ---- Number of potential interactions ---- $
    celltypes_oi <- c(all_df %>%
        pull(id) %>%
        unique())
    celltype_pairs <- tidyr::crossing(source = celltypes_oi, target = celltypes_oi) %>% mutate(source_target = paste0(source, "__", target))

    samples_oi <- all_df %>%
        pull(Sample) %>%
        unique()

    n_possible_interactions <- do.call(rbind, lapply(samples_oi, function(sample) {
        # Using all_df = average gene expression
        df_sub <- all_df %>% filter(Sample == sample)
        celltypes_subset <- celltype_pairs %>% filter(source %in% df_sub$id, target %in% df_sub$id)
        celltypes_subset$Sample <- sample
        n_interactions_by_celltype_pair <- sapply(seq_len(nrow(celltypes_subset)), function(i) {
            source_genes <- df_sub %>%
                filter(id == celltypes_subset[i, "source"] %>% pull()) %>%
                pull(expressed_genes) %>%
                str_split(", ", simplify = TRUE) %>%
                unlist() %>%
                unique()
            target_genes <- df_sub %>%
                filter(id == celltypes_subset[i, "target"] %>% pull()) %>%
                pull(expressed_genes) %>%
                str_split(", ", simplify = TRUE) %>%
                unlist() %>%
                unique()

            n_interactions <- interactions_db %>%
                filter(ligand %in% source_genes, receptor %in% target_genes) %>%
                pull(complex_interaction) %>%
                unique() %>%
                length()
            return(n_interactions)
        })
        celltypes_subset$n_possible_interactions <- n_interactions_by_celltype_pair
        return(celltypes_subset)
    }))

    all_n_interactions_by_sample <- all_n_interactions_by_sample %>% left_join(n_possible_interactions %>% select(Sample, source_target, n_possible_interactions))
    saveRDS(all_n_interactions_by_sample, glue("{args$output_dir}/all_n_interactions_by_sample_{type_of_voting}.rds"))

    plt_possible_interactions <- ggplot(data = all_n_interactions_by_sample) +
        geom_point(aes(x = n_possible_interactions, y = n_interactions, color = target)) +
        facet_grid(Region_Grouped ~ source) +
        default_theme() +
        labs(x = "Number of possible interactions", y = "Number of interactions", title = "Number of interactions per cell type pair", subtitle = glue("min. fraction of cells that express the ligand/receptor > {args$min_pct}")) +
        stat_cor(aes(x = n_possible_interactions, y = n_interactions), method = "spearman")
    ggsave(plot = plt_possible_interactions, filename = glue("{args$output_dir}/ninteractions_per_celltypepair_{type_of_voting}_minpct{args$min_pct * 100}.pdf"), width = 25, height = 10, dpi = 300)


    plot_possible_interactions_overall <- ggplot(data = all_n_interactions_by_sample) +
        geom_point(aes(x = n_possible_interactions, y = n_interactions)) +
        default_theme() +
        facet_grid(~Region_Grouped) +
        labs(x = "Number of possible interactions", y = "Number of interactions", title = "Number of interactions per cell type pair", subtitle = glue("min. fraction of cells that express the ligand/receptor > {args$min_pct}")) +
        stat_cor(aes(x = n_possible_interactions, y = n_interactions), method = "spearman")
    ggsave(plot = plot_possible_interactions_overall, filename = glue("{args$output_dir}/ninteractions_per_celltypepair_{type_of_voting}_minpct{args$min_pct * 100}__overall.pdf"), width = 15, height = 10, dpi = 300)
}
