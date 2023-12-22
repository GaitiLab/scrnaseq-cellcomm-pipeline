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
        description = "Post-processing interaction weights from RaCInG",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/RaCInG")
    args$pair <- "tumor_tam"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

log_info("Obtain paths for all patients...")
weights <- list.files(glue("{here::here()}/output/RaCInG/weights__{args$pair}"), full.names = TRUE)
log_info(glue("Number of weights files: {length(weights)}"))

log_info("Format weights...")
all_weights <- do.call(rbind, lapply(weights, function(file) {
    patient_id <- str_split(get_name(file), pattern = "__", simplify = TRUE)[1]
    patient_weights <- read.table(file,
        header = TRUE, sep = ",",
        row.names = 1
    ) %>%
        rownames_to_column("ligand") %>%
        pivot_longer(cols = !ligand, names_to = "receptor", values_to = "weight") %>%
        mutate(patient_id = patient_id)
    return(patient_weights)
}))

saveRDS(all_weights, glue("{args$output_dir}/all_weights__{args$pair}.rds"))
log_info("COMPLETED!")
