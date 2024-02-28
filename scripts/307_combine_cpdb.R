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
        description = "Combine downsampling runs CellPhoneDB per sample"
    )
    parser$add_argument("--input_dir", type = "character", help = "Directory where all runs for the sample can be found", default = "")
    parser$add_argument("--sample_id", type = "character", help = "Identifier for sample", default = "")
    parser$add_argument("--n_repeats", type = "numeric", help = "N times downsampling", default = 2)
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/test_pipeline/303_postproc_cpdb")
    args$input_dir <- glue("{here::here()}/output/test_pipeline/303_postproc_cpdb")
    args$sample_id <- "6419_cortex"
    args$n_repeats <- 3
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
pacman::p_load(harmonicmeanp)

log_info(glue("Load and combine all files for sample: {args$sample_id}"))
obj_list <- list.files(args$input_dir,
    full.names = TRUE,
    pattern = glue("cpdb__{args$sample_id}__run__\\d__postproc\\.rds")
)
all_obj <- do.call(rbind, lapply(obj_list, readRDS))

log_info("Combine p-values (harmonicmeanp) and scores (mean)...")
all_obj <- all_obj %>%
    mutate(
        # p.hmp unable to handle pval of exactly zero, so setting it to the smallest number possible
        #  e.g.     `p.hmp(0)` results in
        #           Error in tailsEstable(x, stableParamObj) :
        #           NA/NaN/Inf in foreign function call (arg 7)
        pval = ifelse(pval == 0, .Machine$double.xmin, pval)
    ) %>%
    group_by(method, Sample, source_target, simple_interaction, complex_interaction) %>%
    reframe(
        # Only combine p-values if number of p-values is equal to number of repeats (number of times of downsampling)
        pval = ifelse(n() == args$n_repeats,
            p.hmp(p = pval, w = NULL, L = n()), NA
        ),
        interaction_score = ifelse(n() == args$n_repeats, mean(interaction_score), NA), n_detected = n()
    )
head(all_obj)
log_info("Save RDS...")
saveRDS(all_obj, glue("{args$output_dir}/cpdb__{args$sample_id}__postproc.rds"))
log_info("COMPLETED!")
