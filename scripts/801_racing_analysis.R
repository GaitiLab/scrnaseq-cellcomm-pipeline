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
        description = "Analyse RaCInG results",
    )
    parser$add_argument("-p", "--pair", help = "Pair of cell types", type = "character")
    parser$add_argument("-t", "--top_n", help = "Number of top interactions to select", type = "integer")
    parser$add_argument("-m", "--min_frac", help = "Minimum fraction of patients with interaction", type = "numeric", default = 0)
    parser$add_argument("-n", "--min_patients", help = "Minimum number of patients with interaction", type = "integer", default = 0)
    parser$add_argument("-s", "--n_sd", help = "Number of SDs above mean to select interactions", type = "numeric", default = 0)
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/RaCInG/figures")
    args$pair <- "tam_tumor"

    args$top_n <- 25
    args$min_frac <- 1
    args$min_patients <- 0
    args$n_sd <- 0
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
pacman::p_load_gh("jokergoo/ComplexHeatmap", "jmw86069/colorjam")
pacman::p_load(circlize, ggplot2, ggrepel, ggpubr, rstatix, VennDiagram)

log_info("Load metadata...")
metadata <- readRDS(glue("{here::here()}/output/RaCInG/GBM_Neftel2019_subtypes.rds")) %>% mutate(patient_id = substr(str_replace_all(sample_id, "\\.", "-"), 1, 15))
head(metadata)

log_info("Load objects with interaction weights...")
all_weights <- readRDS(glue("{here::here()}/output/RaCInG/all_weights__{args$pair}.rds"))

log_info("Add metadata to weights...")
all_weights <- all_weights %>%
    left_join(metadata, by = "patient_id") %>%
    mutate(interaction = paste(ligand, receptor, sep = " - "))
head(all_weights)

# Get unique subtypes
subtypes <- all_weights %>%
    pull(Neftel3) %>%
    unique()

# Select the top interactions for each subtype
top_interactions <- do.call(rbind, lapply(subtypes, function(subtype) {
    n_patients <- all_weights %>%
        filter(Neftel3 == subtype) %>%
        pull(patient_id) %>%
        unique() %>%
        length()
    log_info(glue("Patients with {subtype}: {n_patients}"))

    mean_weights <- all_weights %>%
        filter(Neftel3 == subtype) %>%
        mutate(is_present = ifelse(weight > 0, 1, 0)) %>%
        group_by(interaction) %>%
        summarise(mean = mean(weight), total_patients = n_patients, n_present = sum(is_present), frac = sum(n_present) / total_patients) %>%
        ungroup() %>%
        arrange(desc(mean))

    weight_mean <- mean_weights %>%
        pull(mean) %>%
        mean()
    weight_sd <- mean_weights %>%
        pull(mean) %>%
        sd()

    log_info(glue("Mean weight: {weight_mean}"))
    log_info(glue("SD weight: {weight_sd}"))

    mean_weights <- mean_weights %>%
        filter(
            # mean > weight_mean + args$n_sd * weight_sd,
            frac >= args$min_frac,
            n_present >= args$min_patients
        ) %>%
        arrange(desc(mean)) %>%
        head(args$top_n) %>%
        mutate(top_interactions = subtype)
    return(mean_weights)
}))
unique_interactions <- top_interactions %>%
    pull(interaction) %>%
    unique()
log_info(glue("Number of unique interactions: {length(unique_interactions)}"))

saveRDS(top_interactions, glue("{here::here()}/output/RaCInG/top_interactions__{args$pair}.rds"))
# hm_df <- all_weights %>%
#     filter(interaction %in% unique_interactions) %>%
#     pivot_wider(id_cols = c(patient_id, Neftel3), names_from = interaction, values_from = weight)

# hm_mat <- hm_df %>% column_to_rownames("patient_id") %>% select(-Neftel3) %>% as.matrix()

# rowAnnotation(subtype = df %>% select(Neftel3) %>% as.matrix())

# create_hm(mat = hm_mat, output_file = glue("{args$output_dir}/hm.pdf"),
# save_plot = TRUE, col_fun = colorRamp2(c(0, mean(hm_mat)), c("white", "red")),
# left_annotation = rowAnnotation(subtype = df %>% select(Neftel3) %>%
# as.matrix()))
log_info("Compare weights between subtypes...")
all_weights_subset <- all_weights %>% filter(interaction %in% unique_interactions)

stat_test <- all_weights_subset %>%
    group_by(interaction) %>%
    wilcox_test(weight ~ Neftel3) %>%
    add_y_position(fun = "median_iqr", step.increase = 0.02)

ggboxplot(all_weights_subset, x = "Neftel3", y = "weight", fill = "Neftel3", facet.by = "interaction", scales = "free_y") + stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, ncol = 3) + rotate_x_text()

ggsave(glue("{args$output_dir}/{args$pair}__boxplot.pdf"), width = 15, height = 20)

log_info("Create Venn diagram...")
venn_data <- top_interactions %>%
    select(-mean, -total_patients, -frac) %>%
    reshape2::dcast(interaction ~ top_interactions, value.var = "interaction") %>%
    column_to_rownames("interaction")

venn_data_ls <- lapply(venn_data %>% as.list(), function(x) x[!is.na(x)])

venn.diagram(
    x = venn_data_ls, category.names = names(venn_data_ls), filename = glue("{args$output_dir}/{args$pair}__venn.png"),
    output = FALSE,
    # main = str_replace_all(pair, "__", " - "),
    disable.logging = TRUE,
    # Output features
    imagetype = "png",
    height = 480,
    width = 500,
    resolution = 600,
    compression = "lzw",

    # Circles
    lwd = 2,
    lty = "blank",
    fill = rainbowJam(length(names(venn_data_ls))),

    # Numbers
    cex = .25,
    fontface = "bold",
    fontfamily = "sans",

    # Set names
    cat.cex = 0.17,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans",

    # Title
    main.fontfamily = "sans",
    main.cex = 0.3,
    main.fontface = "bold",
    margin = 0.05
)

# Get unique interactions
ac_like_unique <- setdiff(venn_data %>% pull("AC-like"), c(venn_data %>% pull("NPC & OPC-like"), venn_data %>% pull("MES-like")))
npc_opc_like_unique <- setdiff(venn_data %>% pull("NPC & OPC-like"), c(venn_data %>% pull("AC-like"), venn_data %>% pull("MES-like")))
mes_like_unique <- setdiff(venn_data %>% pull("MES-like"), c(venn_data %>% pull("AC-like"), venn_data %>% pull("NPC & OPC-like")))

log_info(glue("Unique AC-like interactions: {paste(ac_like_unique, collapse = ', ')}"))
log_info(glue("Unique NPC & OPC-like interactions: {paste(npc_opc_like_unique, collapse = ', ')}"))
log_info(glue("Unique MES-like interactions: {paste(mes_like_unique, collapse = ', ')}"))
