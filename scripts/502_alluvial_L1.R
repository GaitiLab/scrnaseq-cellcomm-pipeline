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
        description = "Create Alluvial/Sankey Diagram"
    )
    parser$add_argument("--input_file", type = "character", help = "Input file")
    parser$add_argument("--interactions_db",
        type = "character", help = "Interactions database (default = data/interactions_db/ref_db.rds)", default =
            "{here::here()/data/interactions_db/ref_db.rds}"
    )
    parser$add_argument("--agg_level", type = "character", help = "either 1) sample or 2) patient, dependent on whether scripts 402a and 402b have been executed (default = Patient)", default = "Patient")
    parser$add_argument("--condition_varname", type = "character", help = "Variable for grouping", default = "Region_Grouped")
    parser$add_argument("--colors", type = "character", help = "Path to RDS file containing a named list for the colors.")
    parser$add_argument("--annot", type = "character", help = "Annotation variable name", default = "CCI_CellClass_L1")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L1"
    run_name <- "CCI_CellClass_L1"
    args$agg_level <- "patient"
    args$output_dir <- glue("{here::here()}/output/{run_name}/502_alluvial")
    args$input_file <- glue("{here::here()}/output/{run_name}/402_aggregation/402_{args$agg_level}_interactions_mvoted_w_filters.rds")
    args$interactions_db <- "001_data_local/interactions_db_v2/ref_db.rds"
    args$colors <- glue("{here::here()}/000_misc_local/{args$annot}_network_colors.rds")
    args$condition_varname <- "Region_Grouped"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
output_dir <- glue("{args$output_dir}/{args$agg_level}")
create_dir(output_dir)

# Load additional libraries
log_info("Load additional libraries...")
pacman::p_load(ggplot2, ggalluvial, ggtext)

log_info("Load interactions...")
interactions <- readRDS(args$input_file) %>% mutate(Region_Grouped = factor(Region_Grouped, levels = c("PT", "TE", "SC")))


interactions <- interactions %>%
    mutate(source_target_abbrev = str_replace_all(source_target, CELLTYPE_ANNOT_ABBREV_DICT)) %>%
    separate(source_target, c("source", "target"), sep = "__", remove = FALSE) %>%
    separate(source_target_abbrev, c("source_abbrev", "target_abbrev"), sep = "__")


# Defaults
malignant_name <- "Malignant"
alluv_width <- .15
alluv_linewidth <- .75
alluv_linecolor <- "black"
add_on_theme <- theme(
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), plot.background = element_blank(),
    axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank()
)

palettes <- list(palette1 = COLORS_L1_TME, palette2 = COLORS_L1_TME_HIGHLIGHTED)


for (option in INTERACTIONS_POST_FILTERING_OPTIONS) {
    if (!option %in% colnames(interactions)) {
        next
    }

    mal_to_tme <- interactions %>%
        rowwise() %>%
        filter(
            str_detect(source, malignant_name),
            !(str_detect(source, malignant_name) && str_detect(target, malignant_name)),
            !!sym(option), !is.na(!!sym(option)),
        ) %>%
        group_by(
            !!sym(args$condition_varname),
            source_abbrev, source, target_abbrev, target
        ) %>%
        summarise(Freq = n()) %>%
        as.data.frame() %>%
        rowwise() %>%
        mutate(
            pair_ordered = paste0(sort(c(source_abbrev, target_abbrev)), collapse = "__"),
            source_target = factor(paste0(c(source_abbrev, target_abbrev), collapse = "__"))
        ) %>%
        ungroup()
    #     # A tibble: 6 × 7
    #   Region_Grouped source_abbrev source    target_abbrev target     Freq source_target
    #   <fct>          <chr>         <chr>     <chr>         <chr>     <int> <fct>
    # 1 PT             Astro         Astrocyte Mal           Malignant     8 Astro__Mal
    # 2 PT             MG            Microglia Mal           Malignant    36 MG__Mal
    # 3 PT             Mal           Malignant Astro         Astrocyte     5 Mal__Astro
    # 4 PT             Mal           Malignant MG            Microglia    30 Mal__MG
    # 5 PT             Mal           Malignant N             Neuron       50 Mal__N
    # 6 PT             Mal           Malignant OPC           OPC          14 Mal__OPC

    tme_to_mal <- interactions %>%
        rowwise() %>%
        filter(str_detect(target, malignant_name), !!sym(option), !(str_detect(source, malignant_name) && str_detect(target, malignant_name))) %>%
        group_by(Region_Grouped, source, source_abbrev, target_abbrev) %>%
        summarise(Freq = n()) %>%
        as.data.frame() %>%
        rowwise() %>%
        mutate(
            # pair_ordered = paste0(sort(c(source_abbrev, target_abbrev)), collapse = "__"),
            source_target = factor(paste0(c(source_abbrev, target_abbrev), collapse = "__"))
        ) %>%
        ungroup()
    # A tibble: 6 × 6
    #   Region_Grouped source          source_abbrev target_abbrev  Freq source_target
    #   <fct>          <chr>           <chr>         <chr>         <int> <fct>
    # 1 PT             Astrocyte       Astro         Mal               8 Astro__Mal
    # 2 PT             Microglia       MG            Mal              36 MG__Mal
    # 3 PT             Neuron          N             Mal              48 N__Mal
    # 4 PT             OPC             OPC           Mal              21 OPC__Mal
    # 5 PT             Oligodendrocyte Oligo         Mal              40 Oligo__Mal
    # 6 TE             Macrophage      M             Mal              38 M__Mal
    for (i in names(palettes)) {
        color_palette <- palettes[[i]]
        # Visualizing
        p1 <- ggplot(
            data = data.frame(mal_to_tme),
            aes(y = Freq, axis1 = Region_Grouped, axis2 = source_abbrev, axis3 = target_abbrev)
        ) +
            geom_alluvium(aes(fill = target), width = alluv_width, alpha = 1, linewidth = alluv_linewidth) +
            geom_stratum(alpha = 0.15, width = alluv_width, linewidth = alluv_linewidth) +
            geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
            scale_x_discrete(limits = c("Region", "Tumor\n(Sender)", "TME\n(Receiver)"), expand = c(.05, .05)) +
            scale_fill_manual(values = color_palette) +
            custom_theme() +
            add_on_theme +
            geom_flow(
                stat = "alluvium", aes(fill = target), width = alluv_width,
                color = alluv_linecolor, linewidth = alluv_linewidth
            ) +
            labs(title = "Tumor - TME")

        p2 <- ggplot(
            tme_to_mal,
            aes(y = Freq, axis1 = Region_Grouped, axis2 = source_abbrev, axis3 = target_abbrev)
        ) +
            geom_alluvium(aes(fill = source), width = alluv_width, alpha = 1, linewidth = alluv_linewidth) +
            geom_stratum(alpha = 0.15, width = alluv_width, linewidth = alluv_linewidth) +
            geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
            scale_x_discrete(limits = c("Region", "TME\n(Sender)", "Tumor\n(Receiver)"), expand = c(.05, .05)) +
            scale_fill_manual(values = color_palette) +
            custom_theme() +
            add_on_theme +
            geom_flow(
                stat = "alluvium", aes(fill = source), width = alluv_width,
                color = alluv_linecolor, linewidth = alluv_linewidth
            ) +
            labs(title = "TME - Tumor")


        p <- ggpubr::ggarrange(p1, p2)
        p

        ggsave(glue("{output_dir}/alluvial__{option}__{i}.pdf"), p, width = 15, height = 8)
    }
}
