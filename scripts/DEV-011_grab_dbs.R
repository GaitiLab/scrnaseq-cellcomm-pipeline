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
    args$output_dir <- glue("{here::here()}/001_data_local/interactions_db_v2/")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

old_liana <- readRDS(glue("{here::here()}/001_data_local/interactions_db/custom_liana.rds"))

# Load additional libraries
# pacman::p_load("jinworks/CellChat")
pacman::p_load("CellChat")
pacman::p_load_gh("Wei-BioMath/NeuronChat")
# pacman::p_load("liana")

# Base Database (LIANA)
log_info("Loading LIANA database: Consensus...")
liana_consensus <- liana::select_resource("Consensus")[[1]] %>% mutate(method = "LIANA consensus")


log_info("Loading LIANA database: Ramilowski2015...")
liana_ramilowski <- liana::select_resource("Ramilowski2015")[[1]] %>% mutate(method = "LIANA Ramilowski2015")

log_info("Account for missing columns...")
missing_liana <- setdiff(colnames(liana_consensus), colnames(liana_ramilowski))
missing_ramilowski <- setdiff(colnames(liana_ramilowski), colnames(liana_consensus))

liana_consensus[missing_ramilowski] <- NA
liana_ramilowski[missing_liana] <- NA

log_info("Combining the two LIANA databases...")
liana_db <- rbind(liana_consensus, liana_ramilowski) %>%
    liana::decomplexify() %>%
    mutate(interaction = paste0(source_genesymbol_complex, "__", target_genesymbol_complex))


log_info(glue("Number of interactions: {nrow(liana_consensus)}"))
log_info(glue("Number of interactions: {nrow(liana_ramilowski)}"))
log_info(glue("Number of interactions: {nrow(liana_db %>% distinct(interaction))}"))

#----- CellChat V2 ----- #
cellchat_db <- CellChatDB.human
cellchat_db_interaction <- cellchat_db$interaction
cellchat_db_complex <- cellchat_db$complex
cellchat_db_cofactor <- cellchat_db$cofactor
cellchat_db_geneinfo <- cellchat_db$geneInfo

cellchat_custom <- cellchat_db$interaction %>%
    select(ligand.symbol, receptor.symbol, interaction_name_2, ligand, receptor) %>%
    mutate(
        source_genesymbol = ifelse(ligand.symbol == "", ligand, str_split(ligand.symbol, ",", simplify = TRUE)[, 1]),
        target_genesymbol = ifelse(receptor.symbol == "", receptor, str_split(receptor.symbol, ",", simplify = TRUE)[, 1]),
        source_genesymbol_complex = ifelse(ligand.symbol == "", str_replace_all(ligand, ", ", "_"), str_replace_all(ligand.symbol, ", ", "_")),
        target_genesymbol_complex = ifelse(receptor.symbol == "", str_replace_all(receptor, ", ", "_"), str_replace_all(receptor.symbol, ", ", "_")),
        method = "CellChatDB v2",
    ) %>%
    select(-ligand.symbol, -receptor.symbol, -interaction_name_2, -ligand, -receptor) %>%
    mutate(interaction = paste0(source_genesymbol_complex, "__", target_genesymbol_complex)) %>%
    left_join(cellchat_db_geneinfo %>% select("EntryID.uniprot", "Symbol"), by = c("source_genesymbol" = "Symbol")) %>%
    rename(source = "EntryID.uniprot") %>%
    left_join(cellchat_db_geneinfo %>% select("EntryID.uniprot", "Symbol"), by = c("target_genesymbol" = "Symbol")) %>%
    rename(target = "EntryID.uniprot")

# cellchat_db_complex[cellchat_db_complex == ""] <- NA
# cellchat_db_complex <- cellchat_db_complex %>% unite(., col = "genes_complex", na.rm = TRUE, sep = "_")

# out <- cellchat_custom %>% left_join(cellchat_db_complex %>% rownames_to_column("complex"), by = c("receptor" = "complex")) %>% select(target_genesymbol_complex, genes_complex)

#  na.rm = TRUE, sep = ", ", remove = FALSE
rownames(cellchat_custom) <- NULL
head(cellchat_custom)

# CD94:NKG2A

missing_cellchat <- setdiff(colnames(liana_db), colnames(cellchat_custom))
cellchat_custom[missing_cellchat] <- NA

log_info("Add CellChat to LIANA database...")
liana_db <- rbind(liana_db, cellchat_custom)

log_info(glue("Number of interactions: {nrow(cellchat_custom)}"))
log_info(glue("Number of interactions: {nrow(liana_db %>% distinct(interaction))}"))
# source, target, source_genesymbol, source_genesymbol_complex,
# target_genesymbol_complex, target_genesymbol

# CellPhoneDB V5
cpdb_complex_input <- read.csv(glue("{here::here()}/001_data_local/interactions_db_v2/cellphonedb_v5.0.0/complex_input.csv"), sep = ",")

cpdb_gene_input <- read.csv(glue("{here::here()}/001_data_local/interactions_db_v2/cellphonedb_v5.0.0/gene_input.csv"), sep = ",")

cpdb_interaction_input <- read.csv(glue("{here::here()}/001_data_local/interactions_db_v2/cellphonedb_v5.0.0/interaction_input.csv"), sep = ",")

cpdb_protein_input <- read.csv(glue("{here::here()}/001_data_local/interactions_db_v2/cellphonedb_v5.0.0/protein_input.csv"), sep = ",")

custom_cellphonedb <- cpdb_interaction_input %>%
    select(partner_a, partner_b, interactors) %>%
    rename(source = partner_a, target = partner_b)

interaction_units <- str_split(custom_cellphonedb$interactors, "-", simplify = TRUE)
colnames(interaction_units) <- paste0("unit_", seq_len(ncol(interaction_units)))

custom_cellphonedb <- cbind(custom_cellphonedb, interaction_units) %>%
    mutate(
        source_genesymbol_complex = case_when(nchar(unit_2) == 1 ~ paste0(unit_1, "-", unit_2), TRUE ~ unit_1),
        target_genesymbol_complex = case_when(
            nchar(unit_2) == 1 ~ unit_3,
            nchar(unit_3) == 1 ~ paste0(unit_2, "-", unit_3), TRUE ~ unit_2
        ), method = "CellPhoneDB v5"
    ) %>%
    select(-all_of(colnames(interaction_units)))


ligand <- str_split(custom_cellphonedb$source_genesymbol_complex, "\\+", simplify = TRUE)[, 1]
receptor <- str_split(custom_cellphonedb$target_genesymbol_complex, "\\+", simplify = TRUE)[, 1]

custom_cellphonedb$source_genesymbol <- ligand
custom_cellphonedb$target_genesymbol <- receptor

custom_cellphonedb <- custom_cellphonedb %>%
    select(source, target, source_genesymbol, target_genesymbol, source_genesymbol_complex, target_genesymbol_complex, method) %>%
    mutate(source_genesymbol_complex = str_replace_all(source_genesymbol_complex, "\\+", "_"), target_genesymbol_complex = str_replace_all(target_genesymbol_complex, "\\+", "_"))
saveRDS(custom_cellphonedb, glue("{args$output_dir}/TMP_cellphonedb__{get_current_date()}.rds"))
missing_cellphonedb <- setdiff(colnames(liana_db), colnames(custom_cellphonedb))

custom_cellphonedb[missing_cellphonedb] <- NA

log_info("Add CellPhoneDB to LIANA database...")
liana_db <- rbind(liana_db, custom_cellphonedb)

log_info(glue("Number of interactions: {nrow(custom_cellphonedb)}"))
log_info(glue("Number of interactions: {nrow(liana_db %>% distinct(interaction))}"))

liana_db_filtered <- liana_db %>%
    filter(!is.na(interaction)) %>%
    distinct(pick("interaction", "source_genesymbol_complex", "target_genesymbol_complex"), .keep_all = TRUE) %>%
    mutate(simple_interaction = paste0(source_genesymbol, "__", target_genesymbol))

to_remove <- liana_db_filtered %>%
    filter(!(simple_interaction == interaction)) %>%
    pull(simple_interaction) %>%
    unique()

liana_db_final <- liana_db_filtered %>% filter(!((simple_interaction %in% to_remove) & (simple_interaction == interaction)))
log_info(glue("Number of interactions: {nrow(liana_db_final %>% distinct(interaction))}"))



saveRDS(liana_db_final, glue("{args$output_dir}/liana_db__{get_current_date()}.rds"))
write.csv(liana_db_final, glue("{args$output_dir}/liana_db__{get_current_date()}.csv"), row.names = FALSE)
log_info("COMPLETED!")
