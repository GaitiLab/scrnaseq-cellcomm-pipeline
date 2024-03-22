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
        description = "Get metadata",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    run_name <- "CCI_CellClass_L4_w_agg"
    args$annot <- "CCI_CellClass_L2_2"

    args$output_dir <- glue("{here::here()}/output/{run_name}/504_hockey")
    args$input_file <- glue("{here::here()}/output/{run_name}/402_aggregation/402_interactions_combi_agg.rds")
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
pacman::p_load(ggplot2, ggrepel)
malignant_is_sender <- TRUE

# ognize malignant cells
if (args$annot == "CCI_CellClass_L1") {
    malign_id <- "Malignant"
} else {
    malign_id <- "_like|-like"
}
# malign_function <- ifelse(malignant_is_sender, "source", "target")

obj <- readRDS(args$input_file) %>%
    separate(source_target, c("source", "target"), sep = "__", remove = FALSE) %>%
    filter(
        # Only interactions between Malignant - TME
        str_detect(source_target, malign_id),
        !str_count(source_target, malign_id) == 2,

        # No malignant-malignant interactions
        # !(str_detect(source, malign_id) && str_detect(target, malign_id))
    ) %>%
    mutate(complex_interaction = str_replace_all(complex_interaction, "__", "-"), source_target = str_replace_all(source_target, "__", "\n"))

head(obj)

cell_type_pairs <- obj %>%
    pull(source_target) %>%
    unique()


tmp <- obj %>%
    mutate(log10p = -log10(pval)) %>%
    group_by(Region_Grouped, source_target, .add = TRUE) %>%
    # arrange(desc(log10p), .by_group = TRUE) %>%
    mutate(order_id = min_rank(desc(log10p)))


ggplot(data = tmp %>% filter(source_target %in% cell_type_pairs[1:5]), aes(x = order_id, y = log10p, color = Region_Grouped)) +
    geom_point() +
    custom_theme() +
    labs(y = "-log10(pval)", x = "rank by -log10(pval)") +
    geom_label_repel(aes(label = ifelse(order_id < 5, complex_interaction, "")), min.segment.length = unit(0.1, "lines")) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text.y.right = element_text(angle = 0, hjust = 0)) +
    facet_grid(source_target ~ Region_Grouped, scales = "free_x")

ggsave(filename = glue("{args$output_dir}/hockey_plot.pdf"), width = 14, height = 45)
