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
        description = "Update CellPhoneDB database w/ LIANA",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/001_data_local/interactions_db_v2/cellphonedb_custom/")
    args$liana_db <- "001_data_local/interactions_db_v2/liana_db_updated_final.rds"
    args$cpdb <- glue("{here::here()}/001_data_local/interactions_db_v2/cellphonedb_v5.0.0/")
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
# mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") # human
# gene_ref <- biomaRt::getBM(
#     attributes = c("ensembl_gene_id", "uniprotswissprot", "hgnc_symbol", "entrezgene_id"),
#     mart = mart
# )
# saveRDS(gene_ref, "001_data_local/interactions_db_v2/biomart_hsapiens.rds")
gene_ref <- readRDS("001_data_local/interactions_db_v2/biomart_hsapiens.rds")
# up <- UniProt.ws::UniProt.ws(taxId = 9606)
# keys <- gene_ref %>% pull(entrezgene_id) %>% unique()
# columns <- c("accession", "id")
# kt <- "GeneID"
# res <- UniProt.ws::select(up, keys, columns, kt)
# saveRDS(res, "001_data_local/interactions_db_v2/uniprot_entrez_biomart.rds")
uniprot_proteins <- readRDS("001_data_local/interactions_db_v2/uniprot_entrez_biomart.rds") %>%
    rename(uniprot = Entry, protein_name = Entry.Name, entrezgene_id = From)

gene_protein_ref <- merge(gene_ref, uniprot_proteins, by = "entrezgene_id") %>%
    select(-uniprotswissprot) %>%
    distinct() %>%
    rename(ensembl = ensembl_gene_id) %>%
    distinct(uniprot, .keep_all = TRUE) %>%
    distinct(hgnc_symbol, .keep_all = TRUE)
log_info(glue("Gene-protein reference table: {nrow(gene_protein_ref)}"))

# Manual curation
gene_protein_ref <- gene_protein_ref %>% filter(!(uniprot %in% c("P46091", "P16619")))

log_info("Load LIANA database...")
liana_db <- readRDS(args$liana_db)

log_info("Load CellPhoneDB database files...")
complex_input <- read.csv(glue("{args$cpdb}/complex_input.csv"))
gene_input <- read.csv(glue("{args$cpdb}/gene_input.csv"))
interaction_input <- read.csv(glue("{args$cpdb}/interaction_input.csv"))
protein_input <- read.csv(glue("{args$cpdb}/protein_input.csv"))

interaction_units <- str_split(interaction_input$interactors, "-", simplify = TRUE)
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
    select(source_genesymbol, target_genesymbol) %>%
    mutate(complex_interaction = paste0(source_genesymbol, "_", target_genesymbol))

log_info("Remove interactions that are already present in CellPhoneDB...")
liana_db_filtered <- liana_db %>% filter(!(complex_interaction %in% custom_cellphonedb$complex_interaction))

# Gene input
# mandatory: gene_name, uniprot, hgnc_symbol and ensembl
log_info("Gene input...")
liana_genes <- do.call(c, liana_db_filtered %>%
    select(starts_with("ligand_subunit") | starts_with("receptor_subunit")) %>%
    c())
liana_genes <- liana_genes[(!is.na(liana_genes)) & (liana_genes != "")]

genes_not_in_ref <- setdiff(liana_genes, gene_protein_ref$hgnc_symbol)
liana_genes_filtered <- intersect(liana_genes, gene_protein_ref$hgnc_symbol)
interactions_to_keep <- rowSums(apply(liana_db_filtered %>%
    select(starts_with("ligand_subunit") | starts_with("receptor_subunit")), 2, function(col) {
    col %in% c(genes_not_in_ref)
})) == 0
liana_db_filtered <- liana_db_filtered[interactions_to_keep, ]
log_info(glue("LIANA database: {nrow(liana_db_filtered)}"))

genes_to_add <- setdiff(liana_genes_filtered, gene_input$gene_name)
genes_not_in_ref_tbl <- setdiff(genes_to_add, gene_protein_ref$hgnc_symbol)
log_info(glue("Genes not in reference table: {length(genes_not_in_ref_tbl)}"))
log_info(glue("Genes in reference table: {length(intersect(genes_to_add, gene_protein_ref$hgnc_symbol))}"))

genes_to_add <- data.frame(gene_protein_ref) %>%
    filter(hgnc_symbol %in% genes_to_add) %>%
    mutate(gene_name = hgnc_symbol) %>%
    select(gene_name, uniprot, hgnc_symbol, ensembl)
log_info(glue("Genes to add: {nrow(genes_to_add)}"))
genes_to_add <- genes_to_add %>% distinct(uniprot, .keep_all = TRUE)
genes_to_add <- genes_to_add[colnames(gene_input)]
log_info(glue("Genes to add after removing duplicates: {nrow(genes_to_add)}"))

gene_input_updated <- rbind(gene_input, genes_to_add) %>% distinct()
log_info(glue("Gene input: {nrow(gene_input)}"))
log_info(glue("Gene input: {nrow(gene_input_updated)}"))

# Protein input
# mandatory: uniprot, protein_name
proteins_to_add <- data.frame(gene_protein_ref) %>%
    filter(hgnc_symbol %in% genes_to_add$hgnc_symbol) %>%
    mutate(gene_name = hgnc_symbol) %>%
    select(uniprot, protein_name)

missing_cols <- setdiff(colnames(protein_input), colnames(proteins_to_add))
proteins_to_add[missing_cols] <- ""
proteins_to_add <- proteins_to_add[colnames(protein_input)]
protein_input_updated <- rbind(protein_input, proteins_to_add)
log_info(glue("Protein input before: {nrow(protein_input)}"))
log_info(glue("Protein input after update: {nrow(protein_input_updated)}"))
protein_input_updated <- protein_input_updated %>% distinct(uniprot, .keep_all = TRUE)
log_info(glue("Protein input: {nrow(protein_input_updated)}"))

# Complex input
# mandatory: complex_name, complex_subunits
log_info("Complex input...")
ligand_subunits <- liana_db_filtered %>%
    filter(ligand_subunit_2 != "", !is.na(ligand_subunit_2)) %>%
    select(starts_with("ligand_subunit"))
ligand_protein_subunits <- data.frame(do.call(cbind, lapply(colnames(ligand_subunits), function(col) {
    protein_names <- ligand_subunits %>%
        left_join(gene_input_updated, by = setNames("hgnc_symbol", col), multiple = "first") %>%
        pull(uniprot)
    protein_names[is.na(protein_names)] <- ""
    return(protein_names)
})))
colnames(ligand_protein_subunits) <- paste0("uniprot_", 1:5)
ligand_protein_subunits["complex_name"] <- liana_db_filtered %>%
    filter(ligand_subunit_2 != "", !is.na(ligand_subunit_2)) %>%
    pull(source_genesymbol)

receptor_subunits <- liana_db_filtered %>%
    filter(receptor_subunit_2 != "", !is.na(receptor_subunit_2)) %>%
    select(starts_with("receptor_subunit"))
receptor_protein_subunits <- data.frame(do.call(cbind, lapply(colnames(receptor_subunits), function(col) {
    protein_names <- receptor_subunits %>%
        left_join(gene_input_updated, by = setNames("hgnc_symbol", col), multiple = "first") %>%
        pull(uniprot)
    protein_names[is.na(protein_names)] <- ""
    return(protein_names)
})))
colnames(receptor_protein_subunits) <- paste0("uniprot_", 1:5)
receptor_protein_subunits["complex_name"] <- liana_db_filtered %>%
    filter(simple_interaction != complex_interaction, receptor_subunit_2 != "", !is.na(receptor_subunit_2)) %>%
    pull(target_genesymbol)

complexes_to_add <- rbind(ligand_protein_subunits, receptor_protein_subunits)
missing_cols <- setdiff(colnames(complex_input), colnames(complexes_to_add))
complexes_to_add[missing_cols] <- ""
log_info(glue("Complex input to add before: {nrow(complexes_to_add)}"))
complexes_to_add <- complexes_to_add %>% distinct()
complexes_to_add <- complexes_to_add[colnames(complex_input)]
log_info(glue("Complexes to add after removing duplicates: {nrow(complexes_to_add)}"))

complex_input_updated <- rbind(complex_input, complexes_to_add) %>%
    arrange(desc(version), desc(uniprot_5)) %>%
    distinct(complex_name, .keep_all = TRUE) %>%
    distinct(uniprot_1, uniprot_2, uniprot_3, uniprot_4, uniprot_5, .keep_all = TRUE)
log_info(glue("Complex input before: {nrow(complex_input)}"))
log_info(glue("Complexes to add: {nrow(complexes_to_add)}"))
log_info(glue("Complex input after update: {nrow(complex_input_updated)}"))


# Interaction input
# mandatory:  “partner_a”; “partner_b”; “annotation_strategy”; “source”
log_info("Interaction input...")

interactions_to_add <- liana_db_filtered %>%
    select(-source) %>%
    rename(source = sources) %>%
    left_join(gene_input_updated, by = setNames("hgnc_symbol", "ligand_subunit_1"), multiple = "first") %>%
    rename(partner_a = uniprot) %>%
    left_join(gene_input_updated, by = setNames("hgnc_symbol", "receptor_subunit_1"), multiple = "first") %>%
    rename(partner_b = uniprot) %>%
    mutate(
        partner_a = ifelse((ligand_subunit_2 == "") | (is.na(ligand_subunit_2)), partner_a, source_genesymbol),
        partner_b = ifelse((receptor_subunit_2 == "") | (is.na(receptor_subunit_2)), partner_b, target_genesymbol),
        annotation_strategy = "custom"
    ) %>%
    select(partner_a, partner_b, annotation_strategy, source)

log_info(glue("Interactions to add: {nrow(interactions_to_add)}"))
interactions_to_add <- interactions_to_add %>% distinct(partner_a, partner_b)
log_info(glue("Interactions to add after removing duplicates: {nrow(interactions_to_add)}"))

missing_cols <- setdiff(colnames(interaction_input), colnames(interactions_to_add))
interactions_to_add[missing_cols] <- ""
interactions_to_add <- interactions_to_add[colnames(interaction_input)]

interaction_input_updated <- rbind(interaction_input, interactions_to_add)
log_info(glue("Interaction input before: {nrow(interaction_input)}"))
interaction_input_updated <- interaction_input_updated %>% distinct(partner_a, partner_b, .keep_all = TRUE)
log_info(glue("Interaction input after update: {nrow(interaction_input_updated)}"))


# tmp <- data.frame(t(apply(interaction_input_updated %>% select(partner_a, partner_b), 1, function(row) {
#     row[order(row)]
# })))
# to_remove <- tmp[tmp %>% duplicated(), ] %>%
#     mutate(tmp_interaction = paste0(X1, "_", X2)) %>%
#     pull(tmp_interaction)

interaction_input_updated <- interaction_input_updated %>%
    # mutate(tmp_interaction = paste0(partner_a, "_", partner_b)) %>%
    # filter(!(tmp_interaction %in% to_remove)) %>%
    # select(-tmp_interaction) %>%
    distinct() %>%
    filter(partner_a != "", partner_b != "")

# Save updated files
gene_input_updated[is.na(gene_input_updated)] <- ""
protein_input_updated[is.na(protein_input_updated)] <- ""
complex_input_updated[is.na(complex_input_updated)] <- ""
interaction_input_updated[is.na(interaction_input_updated)] <- ""

# interaction_input_updated %>% distinct(partner_a, partner_b)

# included_proteins <- interaction_input_updated %>%
#     select(partner_a, partner_b) %>%
#     flatten() %>%
#     unlist() %>%
#     unique()

# protein_input_updated <- protein_input_updated %>% filter(uniprot %in% included_proteins)


complex_input_updated <- complex_input_updated %>% mutate(version = ifelse(version == "", "CellPhoneDBcore4.1", version))
log_info("Save updated files...")
write.csv(gene_input_updated, glue("{args$output_dir}/gene_input.csv"), row.names = FALSE)
write.csv(protein_input_updated, glue("{args$output_dir}/protein_input.csv"), row.names = FALSE)

write.csv(complex_input_updated, glue("{args$output_dir}/complex_input.csv"), row.names = FALSE)

write.csv(interaction_input_updated, glue("{args$output_dir}/interaction_input.csv"), row.names = FALSE)


# included_proteins
