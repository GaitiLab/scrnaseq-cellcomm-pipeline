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
        description = "Filter updated LIANA database",
    )
    parser$add_argument("--ref_dir",
        default = "",
        type = "character", help = "Path to reference directory"
    )
    parser$add_argument("--interactions_dir",
        default = "",
        type = "character", help = "Path to interactions directory"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/001_data_local/interactions_db_v2")
    args$interactions_dir <- glue("{here::here()}/001_data_local/interactions_db_v2")
    args$ref_dir <- glue("{here::here()}/000_misc_local/references")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)


log_info("Load additional libraries...")
# mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") # human
# gene_ref <- biomaRt::getBM(
#     attributes = c("ensembl_gene_id", "uniprotswissprot", "hgnc_symbol", "entrezgene_id"),
#     mart = mart
# )
# saveRDS(gene_ref, "001_data_local/interactions_db_v2/biomart_hsapiens.rds")
gene_ref <- readRDS(glue("{args$ref_dir}/biomart_hsapiens.rds"))
# up <- UniProt.ws::UniProt.ws(taxId = 9606)
# keys <- gene_ref %>% pull(entrezgene_id) %>% unique()
# columns <- c("accession", "id")
# kt <- "GeneID"
# res <- UniProt.ws::select(up, keys, columns, kt)
# saveRDS(res, "001_data_local/interactions_db_v2/uniprot_entrez_biomart.rds")
uniprot_proteins <- readRDS(glue("{args$ref_dir}/uniprot_entrez_biomart.rds")) %>%
    rename(uniprot = Entry, protein_name = Entry.Name, entrezgene_id = From)

direct_uniprot <- read.csv(
    glue("{args$ref_dir}/direct_uniprot_download.tsv"),
    sep = "\t", header = TRUE
) %>%
    rename(
        uniprot = Entry,
        protein_name = Entry.Name,
        gene_name = Gene.Names..primary.
    ) %>%
    mutate(ensembl = str_split(Ensembl, ";", simplify = TRUE)[, 1]) %>%
    select(-Ensembl, -Protein.names, -Reviewed)

gene_protein_ref <- merge(gene_ref, uniprot_proteins, by = "entrezgene_id") %>%
    select(-uniprotswissprot) %>%
    distinct() %>%
    rename(ensembl = ensembl_gene_id) %>%
    distinct(uniprot, .keep_all = TRUE) %>%
    distinct(hgnc_symbol, .keep_all = TRUE)
log_info(glue("Gene-protein reference table: {nrow(gene_protein_ref)}"))

cpdb_complex_protein_to_symbol <- readRDS(glue("{args$interactions_dir}/cpdb_complex_protein_to_symbol.rds")) %>% select(complex_gene_name, complex_protein)

# Load additional libraries
protein_complexes <- readRDS(glue("{args$ref_dir}/omnipathr/omnipath_protein_complexes.rds"))
ppi <- readRDS(glue("{args$ref_dir}/omnipathr/omnipath_ppi.rds"))
icn <- readRDS(glue("{args$ref_dir}/omnipathr/omnipath_icn.rds"))

log_info("Loading LIANA (consensus + Ramilowski) combined w/ CellChat and CellPhoneDB...")
liana_db <- readRDS(glue("{args$interactions_dir}/liana_db_combined.rds"))

icn <- icn %>% mutate(across(everything(), as.character))

# Add interaction info from OmniPath
liana_db <- liana_db %>%
    # select(source_genesymbol, target_genesymbol) %>%
    left_join(icn)

# Protein-complexes from OmniPath
components_units <- str_split(protein_complexes$components_genesymbols,
    "_",
    simplify = TRUE
)
protein_complex_comp_ordered <- apply(components_units, 1, function(row) {
    ordered_units <- row[order(row)]
    return(paste0(ordered_units[ordered_units != ""], collapse = "__"))
})

# Removing duplicates (same complex, but different order of subunits)
ligand_subunits <- liana_db %>% select(starts_with("ligand_subunit_"))
ligand_subunits_ordered <- apply(ligand_subunits, 1, function(row) {
    ordered_units <- row[order(row)]
    out <- paste0(ordered_units[(ordered_units != "") & (!is.na(ordered_units))], collapse = "__")
    return(out)
})

receptor_subunits <- liana_db %>% select(starts_with("receptor_subunit_"))
receptor_subunits_ordered <- apply(receptor_subunits, 1, function(row) {
    ordered_units <- row[order(row)]
    out <- paste0(ordered_units[(ordered_units != "") & (!is.na(ordered_units))], collapse = "__")
    return(out)
})

all_ligand_protein_complexes <- protein_complexes
all_ligand_protein_complexes <- all_ligand_protein_complexes %>% mutate(comp_ordered = protein_complex_comp_ordered)
colnames(all_ligand_protein_complexes) <- paste0("lc_", colnames(all_ligand_protein_complexes))
all_receptor_protein_complexes <- protein_complexes
all_receptor_protein_complexes <- all_receptor_protein_complexes %>% mutate(comp_ordered = protein_complex_comp_ordered)
colnames(all_receptor_protein_complexes) <- paste0("rc_", colnames(all_receptor_protein_complexes))

db_tmp <- liana_db %>%
    # select(source_genesymbol, target_genesymbol) %>%
    mutate(ligand_subunits_ord = ligand_subunits_ordered, receptor_subunits_ord = receptor_subunits_ordered) %>%
    left_join(all_ligand_protein_complexes, by = c(ligand_subunits_ord = "lc_comp_ordered")) %>%
    left_join(all_receptor_protein_complexes, by = c(receptor_subunits_ord = "rc_comp_ordered")) %>%
    mutate(n_target_units = str_count(target_genesymbol, "_") + 1, n_ligand_units = str_count(source_genesymbol, "_") + 1)
log_info(glue("Number of interactions: {nrow(db_tmp)}"))
db_tmp <- db_tmp %>% distinct(ligand_subunits_ord, receptor_subunits_ord, .keep_all = TRUE)
log_info(glue("Number of interactions: {nrow(db_tmp)}"))

# Check if there are target units that are not found in reference with PPIs
target_units_not_found <- db_tmp %>% filter(n_target_units > 1, is.na(rc_components))
log_info(glue("Number of target units not found in reference: {nrow(target_units_not_found)}"))
print(head(target_units_not_found))

# Check if there are ligand units that are not found in reference with PPIs
ligand_units_not_found <- db_tmp %>% filter(n_ligand_units > 1, is.na(lc_components))
log_info(glue("Number of ligand units not found in reference: {nrow(ligand_units_not_found)}"))
print(head(ligand_units_not_found))

# Clean up
db_tmp_clean1 <- db_tmp %>% rowwise() %>% 
    mutate(
        source_genesymbol = ifelse(is.na(lc_components_genesymbols) || lc_components_genesymbols == "", source_genesymbol, lc_components_genesymbols),
        target_genesymbol = ifelse(is.na(rc_components_genesymbols) || rc_components_genesymbols == "", target_genesymbol, rc_components_genesymbols),
    ) %>%
    rename(source_protein = lc_components, target_protein = rc_components) %>%
    select(-rc_components_genesymbols, -lc_components_genesymbols, -lc_identifiers, -rc_identifiers, -lc_name, -rc_name) %>% ungroup() 

missing_protein_links <- db_tmp_clean1 %>% filter(
    is.na(source_protein) | is.na(target_protein)
)
log_info(glue("Number of missing protein links: {nrow(missing_protein_links)}"))

# Add proteins that are not found in PPI reference
head(gene_protein_ref)

db_tmp_clean2 <- db_tmp_clean1 %>%
    left_join(gene_protein_ref, by = c(source_genesymbol = "hgnc_symbol")) %>%
    rename(source_tmp_protein_name = uniprot) %>%
    left_join(gene_protein_ref, by = c(target_genesymbol = "hgnc_symbol"), suffix = c("_lc", "_rc")) %>%
    rename(target_tmp_protein_name = uniprot) %>% rowwise() %>% 
    mutate(
        source_protein = ifelse(is.na(source_protein) || source_protein == "", source_tmp_protein_name, source_protein),
        target_protein = ifelse(is.na(target_protein) || target_protein == "", target_tmp_protein_name, target_protein)
    ) %>%
    select(-source_tmp_protein_name, -target_tmp_protein_name) %>% ungroup()

missing_protein_links2 <- db_tmp_clean2 %>% rowwise() %>% filter(
    is.na(source_protein) || is.na(target_protein)
)
log_info(glue("Number of missing protein links: {nrow(missing_protein_links2)}"))

#  Check CPDB complexes
db_tmp_clean3 <- db_tmp_clean2 %>%
    left_join(cpdb_complex_protein_to_symbol, by = c("source_genesymbol" = "complex_gene_name")) %>% rowwise() %>% 
    mutate(source_protein = ifelse(is.na(source_protein) || source_protein == "", complex_protein, source_protein)) %>% ungroup() %>% 
    select(-complex_protein) %>% 
    left_join(cpdb_complex_protein_to_symbol, by = c("target_genesymbol" = "complex_gene_name")) %>% rowwise() %>% 
    mutate(target_protein = ifelse(is.na(target_protein) || target_protein == "", complex_protein, target_protein)) %>%
    select(-complex_protein) %>% ungroup()


missing_protein_links3 <- db_tmp_clean3 %>% rowwise() %>%  filter(
    is.na(source_protein) || is.na(target_protein)
)
log_info(glue("Number of missing protein links: {nrow(missing_protein_links3)}"))
uniprot_proteins <- direct_uniprot %>%
    select(uniprot, gene_name) %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    distinct()

db_tmp_clean4 <- db_tmp_clean3 %>%
    left_join(uniprot_proteins, by = c("source_genesymbol" = "gene_name")) %>% rowwise() %>% 
    mutate(source_protein = ifelse((is.na(source_protein)) || (source_protein == ""), uniprot, source_protein)) %>% ungroup() %>% 
    select(-uniprot) %>%
    left_join(uniprot_proteins, by = c("target_genesymbol" = "gene_name")) %>% rowwise() %>% 
    mutate(target_protein = ifelse(is.na(target_protein) || target_protein == "", uniprot, target_protein)) %>% ungroup() %>% 
    select(-uniprot) %>%
    rowwise() %>%
    mutate(
        source_protein = ifelse(source_protein == "" || is.na(source_protein),
            direct_uniprot[(str_detect(direct_uniprot$Gene.Names, source_genesymbol)), "uniprot"], source_protein
        ),
        target_protein = ifelse(target_protein == "" || is.na(target_protein),
            direct_uniprot[(str_detect(direct_uniprot$Gene.Names, target_genesymbol)), "uniprot"], target_protein
        )
    ) %>%
    ungroup()

missing_links4 <- db_tmp_clean4 %>% filter(
    is.na(source_protein) | is.na(target_protein)
)

log_info(glue("Number of missing protein links: {nrow(missing_links4)}"))
log_info(glue("Number of interactions: {nrow(db_tmp_clean4)}"))

liana_db_updated <- db_tmp_clean4 %>% rowwise() %>% 
    filter(
        !(is.na(source_protein) || is.na(target_protein)), source_protein != "", target_protein != ""
    ) %>% ungroup() %>% 
    distinct(source_genesymbol, target_genesymbol, .keep_all = TRUE)
log_info(glue("Number of interactions: {nrow(liana_db_updated)}"))

ligand_complexes <- liana_db_updated %>%
    select(source_genesymbol, ligand_subunits_ord, source_protein) %>%
    rename(subunits_ordered = ligand_subunits_ord, genesymbol = source_genesymbol, protein = source_protein)
receptor_complexes <- liana_db_updated %>%
    select(target_genesymbol, receptor_subunits_ord, target_protein) %>%
    rename(subunits_ordered = receptor_subunits_ord, genesymbol = target_genesymbol, protein = target_protein)

# Remove duplicates (complexes with same subunits, but different order)
log_info("Removing duplicates...")
complexes <- rbind(ligand_complexes, receptor_complexes)
log_info(glue("Number of complexes: {nrow(complexes)}"))
complexes <- complexes %>%
    distinct(subunits_ordered, protein, .keep_all = TRUE) %>%
    select(-subunits_ordered)
log_info(glue("Number of complexes after removing duplicates: {nrow(complexes)}"))

liana_db_updated <- liana_db_updated %>%
    left_join(complexes, by = c(source_protein = "protein")) %>%
    select(-source_genesymbol) %>%
    rename(source_genesymbol = genesymbol) %>%
    distinct(source_genesymbol, target_genesymbol, .keep_all = TRUE) %>%
    left_join(complexes, by = c(target_protein = "protein")) %>%
    select(-target_genesymbol) %>%
    rename(target_genesymbol = genesymbol) %>%
    distinct(source_genesymbol, target_genesymbol, .keep_all = TRUE)
log_info(glue("Number of interactions: {nrow(liana_db_updated)}"))

log_info("Removing interactions with missing symbols or proteins...")
liana_db_updated <- liana_db_updated %>% rowwise() %>% 
    filter(
        !(source_genesymbol == "" || is.na(source_genesymbol) | target_genesymbol == "" || is.na(target_genesymbol) ||
            source_protein == "" || is.na(source_protein) | target_protein == "" || is.na(target_protein)) 
    ) %>% ungroup()
log_info(glue("Number of interactions: {nrow(liana_db_updated)}"))
# Expected to be the same number of interactions as previous check
liana_db_updated <- liana_db_updated %>%
    separate(source_genesymbol, into = paste0("ligand_subunit_", seq_len(5)), remove = FALSE) %>%
    separate(target_genesymbol, into = paste0("receptor_subunit_", seq_len(5)), remove = FALSE)

# # ------ NEWLY ADDED 24/01/09 -----
# # With distinct first row will be preserved, favoring: (1) LIANA, (2) multi-subunit interaction
# liana_db_updated <- liana_db_updated %>%
#     rowwise() %>%
#     mutate(
#         simple_interaction_ordered = paste0(sort(c(ligand_subunit_1, receptor_subunit_1)), collapse = "__"),
#         is_multi_subunit = (!is.na(ligand_subunit_2) | !is.na(receptor_subunit_2)),
#         method_tmp = case_when(
#             method == "LIANA consensus" ~ 1,
#             method == "LIANA Ramilowski2015" ~ 2,
#             method == "CellphoneDB extracted" ~ 3,
#             method == "CellChat extracted" ~ 3
#         ),
#         n_proof = str_count(references, ";") + str_count(sources, ";")
#     ) %>%
#     arrange(method_tmp, desc(n_proof), desc(is_multi_subunit)) %>%
#     distinct(simple_interaction_ordered, .keep_all = TRUE) %>%
#     ungroup() %>%
#     select(-method_tmp, -simple_interaction_ordered)
# # ------ NEWLY ADDED 24/01/09 -----

ref_db <- liana_db_updated %>% mutate(simple_interaction = paste0(ligand_subunit_1, "__", receptor_subunit_1))

receptor_complex <- apply(liana_db_updated %>% select(starts_with("receptor_subunit_")), 1, function(x) {
    return(paste0(x[(!is.na(x)) & (x != "")], collapse = ":"))
})
ligand_complex <- apply(liana_db_updated %>% select(starts_with("ligand_subunit_")), 1, function(x) {
    return(paste0(x[(!is.na(x)) & (x != "")], collapse = ":"))
})

ref_db <- ref_db %>%
    mutate(
        receptor_complex = receptor_complex,
        ligand_complex = ligand_complex,
        complex_interaction = paste0(ligand_complex, "__", receptor_complex),
        interaction = paste0(source_genesymbol, "_", target_genesymbol)
    ) %>%
    select(source_genesymbol, target_genesymbol, interaction, simple_interaction, complex_interaction, ligand_complex, receptor_complex)
log_info(glue("Number of interactions: {nrow(ref_db)}"))

log_info("Saving reference database...")
# for post-processing/visualization
saveRDS(ref_db, glue("{args$output_dir}/ref_db.rds"))

log_info("Saving LIANA database...")
saveRDS(liana_db_updated, glue("{args$output_dir}/liana_db.rds"))

log_info("Saving LIANA database as CSV for cell2cell...")
write.csv(liana_db_updated, glue("{args$output_dir}/cell2cell_db.csv"), row.names = FALSE)

log_info("COMPLETED!")
