# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
require(GBMutils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Extract cell type signatures from CIBERSORT's signature matrix",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/000_misc_local/")
    args$signatures <- "000_misc_local/genelists.csv"
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
log_info("Load cell type gene signatures...")
signatures <- read.csv(args$signatures, sep = ",")

# Collapse MES (MES1 + MES2) and NPC (NPC1 + NPC2)
Neftel_MES <- c("Neftel_MES1", "Neftel_MES2")
Neftel_NPC <- c("Neftel_NPC1", "Neftel_NPC2")

MES_signature <- signatures %>%
    select(all_of(Neftel_MES)) %>%
    pull()
MES_signature <- MES_signature[MES_signature != ""] %>% unique()

NPC_signature <- signatures %>%
    select(all_of(Neftel_NPC)) %>%
    pull()
NPC_signature <- NPC_signature[NPC_signature != ""] %>% unique()

MES_NPC_signatures_list <- list(Neftel_MES = MES_signature, Neftel_NPC = NPC_signature)

# Deal with AC & OPC
remaining <- c("Neftel_AC", "Neftel_OPC")
remaining_signatures <- signatures %>%
    select(all_of(remaining)) %>%
    as.list()
remaining_signatures_list <- lapply(remaining_signatures, function(list_of_genes) {
    return(list_of_genes[list_of_genes != ""])
})
# Combine
Neftel_signatures <- c(MES_NPC_signatures_list, remaining_signatures_list)

log_info("Save list of signatures...")
saveRDS(Neftel_signatures, file = glue("{args$output_dir}/neftel_signatures.rds"))

log_info("COMPLETED!")
