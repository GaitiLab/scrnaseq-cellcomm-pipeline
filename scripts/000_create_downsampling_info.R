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
        description = "Generate list for downsampling",
    )
    parser$add_argument("--split_varname", type = "character", default = "Sample", help = "Variable for splitting object")
    parser$add_argument("--num_cells", type = "numeric", default = 2500, help = "Number of cells to sample")
    parser$add_argument("--num_repeats", type = "numeric", default = 5, help = "Number of times to sample")
    parser$add_argument("--meta", type = "character", default = "", help = "path to metadata")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
    args$meta <- glue("{here::here()}/output/CCI_CellClass_L1_conf_malign/000_data/gbm_regional_study__metadata.rds")
    args$split_varname <- "Sample"
    args$num_cells <- 2500
    args$num_repeats <- 3
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
meta <- readRDS(args$meta)

# TODO remove later, temporary for testing
samples_oi <- meta %>%
    count(Sample) %>%
    filter(n > 3000) %>%
    pull(Sample) %>%
    unique()
meta <- meta %>% filter(Sample %in% samples_oi)

log_info("Number of samples...")
unique_labels <- meta %>%
    pull(!!sym(args$split_varname)) %>%
    unique()
log_info(glue("N = {length(unique_labels)}"))

sampling_cells <- function(meta, split_varname, label, num_cells) {
    cell_ids <- meta %>%
        filter(!!sym(split_varname) == label) %>%
        rownames()
    cell_ids_subset <- sample(cell_ids, size = num_cells)
    return(cell_ids_subset)
}
log_info("Create run_ids...")
# Create run_ids: {sample_id}__run__{repeat_id}
run_ids <- lapply(unique_labels, function(current_label, num_repeats) {
    paste0(current_label, "__run__", seq_len(num_repeats))
}, num_repeats = args$num_repeats) %>% unlist()

log_info(glue("Sampling {args$num_cells} cells for each sample object and do this {args$num_repeats} times..."))
grid <- pbapply::pblapply(run_ids, function(run_id, metadata, num_cells, split_varname) {
    sample_id <- str_split(run_id, "__", simplify = TRUE)[1]
    return(sampling_cells(meta = metadata, label = sample_id, num_cells = num_cells, split_varname))
}, metadata = meta, num_cells = args$num_cells, split_varname = args$split_varname)
names(grid) <- run_ids

log_info("Save object...")
saveRDS(grid, file = glue("{args$output_dir}/downsampling_info.rds"))

log_info("COMPLETED!")
