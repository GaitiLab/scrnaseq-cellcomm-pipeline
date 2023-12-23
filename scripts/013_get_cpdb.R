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
        description = "Get CellPhoneDB database",
    )
    parser$add_argument("--cpdb",
        type = "character", help = "Path to CellphoneDB database"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/001_data_local/interactions_db_v2/")
    args$cpdb <- "000_misc_local/references/cellphonedb_v5.0.0"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# CellPhoneDB V5
cpdb_interaction_input <- read.csv(glue("{args$cpdb}/interaction_input.csv"), sep = ",")
interaction_units <- str_split(cpdb_interaction_input$interactors, "-", simplify = TRUE)
colnames(interaction_units) <- paste0("unit_", seq_len(ncol(interaction_units)))

custom_cellphonedb <- data.frame(interaction_units) %>%
    mutate(
        source_genesymbol_complex = case_when(nchar(unit_2) == 1 ~ paste0(unit_1, "-", unit_2), TRUE ~ unit_1),
        target_genesymbol_complex = case_when(
            nchar(unit_2) == 1 ~ unit_3,
            nchar(unit_3) == 1 ~ paste0(unit_2, "-", unit_3), TRUE ~ unit_2
        )
    ) %>%
    mutate(source_genesymbol = str_replace_all(source_genesymbol_complex, "\\+", "_"), target_genesymbol = str_replace_all(target_genesymbol_complex, "\\+", "_")) %>%
    select(source_genesymbol, target_genesymbol)
log_info("Saving custom CellPhoneDB...")
saveRDS(custom_cellphonedb, glue("{args$output_dir}/cpdbv5_liana_format.rds"))


log_info("Extract complex output specific for CellPhoneDB...")
complex <- read.csv(glue("{args$cpdb}/complex_input.csv"), sep = ",") %>% mutate(dummy = 1)
db_tmp <- cbind(custom_cellphonedb, cpdb_interaction_input %>% select(partner_a, partner_b))
partners_a <- db_tmp %>%
    select(partner_a, source_genesymbol) %>%
    left_join(complex, by = c("partner_a" = "complex_name")) %>%
    filter(dummy == 1) %>%
    rename(complex_name = partner_a, complex_gene_name = source_genesymbol)
partners_b <- db_tmp %>%
    select(partner_b, target_genesymbol) %>%
    left_join(complex, by = c("partner_b" = "complex_name")) %>%
    filter(dummy == 1) %>%
    rename(complex_name = partner_b, complex_gene_name = target_genesymbol)
cpdb_complex_protein_to_symbol <- rbind(partners_a, partners_b)

cpdb_complex_protein_to_symbol$complex_protein <- apply(cpdb_complex_protein_to_symbol %>% select(starts_with("uniprot_")), 1, function(row) {
    paste0(row[row != "" & !is.na(row)], collapse = "_")
})

log_info("Saving complex output specific for CellPhoneDB...")
saveRDS(cpdb_complex_protein_to_symbol, glue("{args$output_dir}/cpdb_complex_protein_to_symbol.rds"))
log_info("COMPLETED!")
