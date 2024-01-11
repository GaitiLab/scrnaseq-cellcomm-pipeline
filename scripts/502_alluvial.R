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
        description = "Create Alluvial Diagram"
    )
    parser$add_argument("--input_file", type = "character", help = "Input file")
    parser$add_argument("--interactions_db", type = "character", help = "Interactions database")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$annot <- "CCI_CellClass_L2"
    args$output_dir <- glue("{here::here()}/output/{args$annot}/502_alluvial/test_enforce_min2samples")
    args$input_file <- glue("{here::here()}/output/{args$annot}/400_consensus/400_samples_interactions_mvoted_w_filters.rds")
    args$interactions_db <- "001_data_local/interactions_db_v2/ref_db.rds"
    args$meta <- glue("{here::here()}/output/{args$annot}/000_data/gbm_regional_study__metadata.rds")
    args$colors <- glue("{here::here()}/000_misc_local/{args$annot}_network_colors.rds")
    args$has_loops <- FALSE
    args$malignant_name <- "_like"
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
log_info("Load additional libraries...")
pacman::p_load(randomcoloR, ggplot2, ggalluvial)

log_info("Load interactions...")
interactions <- readRDS(args$input_file)

# TODO, comment after generating colors - DONE (DO NOT TOUCH)
# celltypes <- n_cells_per_type[args$annot] %>%
#     unique() %>% pull()
# celltypes <- celltypes[order(celltypes)]
# colors <- distinctColorPalette(length(celltypes))
# names(colors) <- celltypes
# saveRDS(colors, glue("{here::here()}/000_misc_local/{args$annot}_network_colors.rds"))
colors <- readRDS(args$colors)

options <- c("stringent_region", "lenient_region", "lenient_region_pair", "stringent_region_pair")

abbrev_dict <- c(
    "Macrophage" = "M",
    "Microglia" = "MG",
    "Oligodendrocyte" = "Oligo",
    "Astrocyte" = "Astro",
    "Neuron" = "N",
    "OPC" = "OPC",
    "Endothelial" = "Endo",
    "Malignant" = "Mal",
    "Progenitor_like" = "PL",
    "Differentiated_like" = "DL"
)

interactions <- interactions %>%
    mutate(source_target_abbrev = str_replace_all(source_target, abbrev_dict)) %>%
    separate(source_target_abbrev, c("source_abbrev", "target_abbrev"), sep = "__") %>%
    # TODO remove mutate + filter later, just for testing to enforce a min. detection of 2
    mutate(n_detected = str_count(lenient_voting_samples, ",") + 1) %>%
    filter(n_detected >= 2)

for (option in options) {
    # option <- "lenient_region"
    mal_to_other <- interactions %>%
        rowwise() %>%
        filter(str_detect(source, args$malignant_name), !!sym(option), !(str_detect(source, args$malignant_name) && str_detect(target, args$malignant_name))) %>%
        group_by(Region_Grouped, source_abbrev, target_abbrev, target) %>%
        summarise(Freq = n())

    other_to_mal <- interactions %>%
        rowwise() %>%
        filter(str_detect(target, args$malignant_name), !!sym(option), !(str_detect(source, args$malignant_name) && str_detect(target, args$malignant_name))) %>%
        group_by(Region_Grouped, source, source_abbrev, target_abbrev) %>%
        summarise(Freq = n())

    if (args$annot == "CCI_CellClass_L1") {
        p1 <- ggplot(
            as.data.frame(mal_to_other),
            aes(y = Freq, axis2 = target_abbrev, axis1 = Region_Grouped)
        ) +
            geom_alluvium(aes(fill = target), width = 1 / 12) +
            geom_stratum(alpha = .25, width = 1 / 12) +
            geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
            scale_x_discrete(limits = c("Sender", "Receiver"), expand = c(.05, .05)) +
            scale_fill_manual(values = colors) +
            custom_theme() +
            theme(
                legend.position = "none", panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), plot.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
                axis.title.y = element_blank()
            ) +
            labs(title = "Malignant - TME")

        p2 <- ggplot(
            as.data.frame(other_to_mal),
            aes(y = Freq, axis1 = source_abbrev, axis2 = Region_Grouped)
        ) +
            geom_alluvium(aes(fill = source), width = 1 / 12) +
            geom_stratum(alpha = .25, width = 1 / 12) +
            geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
            scale_x_discrete(limits = c("Sender", "Receiver"), expand = c(.05, .05)) +
            scale_fill_manual(values = colors) +
            custom_theme() +
            theme(
                legend.position = "none", panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), plot.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
                axis.title.y = element_blank()
            ) +
            labs(title = "TME - Malignant")
    } else {
        p1 <- ggplot(
            as.data.frame(mal_to_other),
            aes(y = Freq, axis3 = target_abbrev, axis2 = source_abbrev, axis1 = Region_Grouped)
        ) +
            geom_alluvium(aes(fill = target), width = 1 / 12) +
            geom_stratum(alpha = .25, width = 1 / 12) +
            geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
            scale_x_discrete(limits = c("Sender", "Sender - Malignant subtype", "Receiver"), expand = c(.05, .05)) +
            scale_fill_manual(values = colors) +
            custom_theme() +
            theme(
                legend.position = "none", panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), plot.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
                axis.title.y = element_blank()
            ) +
            labs(title = "Malignant - TME")

        p2 <- ggplot(
            as.data.frame(other_to_mal),
            aes(y = Freq, axis1 = source_abbrev, axis2 = target_abbrev, axis3 = Region_Grouped)
        ) +
            geom_alluvium(aes(fill = source), width = 1 / 12) +
            geom_stratum(alpha = .25, width = 1 / 12) +
            geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
            scale_x_discrete(limits = c("Sender", "Sender - Malignant subtype", "Receiver"), expand = c(.05, .05)) +
            scale_fill_manual(values = colors) +
            custom_theme() +
            theme(
                legend.position = "none", panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), plot.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
                axis.title.y = element_blank()
            ) +
            labs(title = "TME - Malignant")
    }




    p <- ggpubr::ggarrange(p1, p2)
    ggsave(glue("{args$output_dir}/alluvial__{option}.pdf"), p, width = 10, height = 5)
    # auto_crop(glue("{args$output_dir}/alluvial__{option}.pdf"))
}
