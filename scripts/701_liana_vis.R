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


obj <- readRDS("~/Desktop/gaitigroup/Users/Joan/nf_work_cci/7c/aa8e2afb1346e54cbafdd4a5dbe6ff/201_cci_liana/liana__6234_2895153_A.rds") %>% liana::liana_aggregate()  %>% filter(aggregate_rank < 0.05) # note that these pvals are already corrected

liana::heat_freq(obj)

obj2 <- readRDS("~/Desktop/gaitigroup/Users/Joan/nf_work_cci/8b/a481da7d26aeb961839a36c41b355e/201_cci_liana/liana__6234_2895153_B.rds") %>% liana::liana_aggregate()  %>% filter(aggregate_rank < 0.05) # note that these pvals are already corrected

liana::heat_freq(obj2)

obj3 <- readRDS("~/Desktop/gaitigroup/Users/Joan/nf_work_cci/94/b0abc8ef93e091807b44ffec3c38dc/201_cci_liana/liana__6467_cortex.rds") %>% liana::liana_aggregate()  %>% filter(aggregate_rank < 0.05) # note that these pvals are already corrected

liana::heat_freq(obj3)


pdf("test.pdf", width = 15, height = 15)
liana::chord_freq(obj3)
dev.off()
# ggplot2::ggsave(filename = glue("{here::here()}/test.png"), width = 10, height = 10, dpi = 300)
