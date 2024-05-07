# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# require(GBMutils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)
if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Update CellPhoneDB database w/ LIANA",
    )
    parser$add_argument("--liana_db",
        default = "",
        type = "character", help = "Path to LIANA database"
    )
    parser$add_argument("--cpdb",
        default = "",
        type = "character", help = "Path to CellPhoneDB database"
    )
    parser$add_argument("--ref_dir",
        default = "",
        type = "character", help = "Path to reference files"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/001_data_local/interactions_db_v2/cellphonedb_custom/")
    args$liana_db <- "001_data_local/interactions_db_v2/liana_db.rds"
    args$cpdb <- glue("{here::here()}/000_misc_local/references/cellphonedb_v5.0.0/")
    args$ref_dir <- glue("{here::here()}/000_misc_local/references/")
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
pacman::p_load("R.utils")

create_dir(glue("{args$output_dir}/sources"))
R.utils::createLink(target = glue("{args$cpdb}/sources"), link = glue("{args$output_dir}/sources"))


log_info("Load CellPhoneDB database files...")
complex_input_cols <- colnames(read.csv(glue("{args$cpdb}/complex_input.csv")))
gene_input_cols <- colnames(read.csv(glue("{args$cpdb}/gene_input.csv")))
interaction_input_cols <- colnames(read.csv(glue("{args$cpdb}/interaction_input.csv")))
protein_input_cols <- colnames(read.csv(glue("{args$cpdb}/protein_input.csv")))

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

direct_uniprot <- read.csv(glue("{args$ref_dir}/direct_uniprot_download.tsv"), sep = "\t", header = TRUE) %>%
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

log_info("Load CellPhoneDB database...")
liana_db <- readRDS(args$liana_db)

# Gene input
# mandatory: gene_name, uniprot, hgnc_symbol and ensemble
source_genes <- liana_db %>%
    separate(source_genesymbol, paste0("ligand_subunit_", seq_len(5)), sep = "_", remove = FALSE) %>%
    select(starts_with("ligand_subunit_")) %>%
    unlist() %>%
    unique()
target_genes <- liana_db %>%
    separate(target_genesymbol, paste0("receptor_subunit_", seq_len(5)), sep = "_", remove = FALSE) %>%
    select(starts_with("receptor_subunit_")) %>%
    unlist() %>%
    unique()
liana_genes <- c(source_genes, target_genes)
liana_genes <- unique(liana_genes[liana_genes != "" & !is.na(liana_genes)])
log_info(glue("Number of genes: {length(liana_genes)}"))

proteins <- str_split(do.call(c, liana_db %>% select(source_protein, target_protein)) %>% unique(), "_", simplify = TRUE)
proteins <- proteins[proteins != ""]
proteins <- data.frame(uniprot = proteins)
proteins <- proteins %>%
    left_join(direct_uniprot, by = "uniprot") %>%
    distinct() %>%
    left_join(gene_protein_ref, by = "uniprot", suffix = c("_x", "_y")) %>%
    mutate(protein_name = ifelse(is.na(protein_name_x) | protein_name_x == "", protein_name_y, protein_name_x)) %>%
    select(-protein_name_x, -protein_name_y) %>%
    rename(ensembl = ensembl_y) %>%
    rowwise() %>%
    mutate(gene_name = str_split(gene_name, ";", simplify = TRUE)[1]) %>%
    mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", hgnc_symbol, gene_name)) %>%
    mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", str_split(Gene.Names, " ", simplify = TRUE)[1], gene_name)) %>%
    ungroup() %>%
    select(gene_name, uniprot, hgnc_symbol, ensembl, protein_name) %>%
    mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "", gene_name, hgnc_symbol)) %>%
    filter(!is.na(gene_name), gene_name != "", hgnc_symbol != "", !is.na(hgnc_symbol), !is.na(ensembl), ensembl != "")
log_info(glue("Number of proteins: {nrow(proteins)}"))

# genes_info <- proteins %>%
#     left_join(uniprot_proteins, by = "uniprot", multiple = "first") %>%
#     distinct()
# genes_info[is.na(genes_info)] <- ""
# log_info(glue("Number of genes: {nrow(genes_info)}"))
missing_genes <- setdiff(liana_genes, proteins$gene_name)
missing_genes2 <- setdiff(proteins$gene_name, liana_genes)
log_info(glue("Number of missing genes: {length(missing_genes)}"))

genes_to_keep <- intersect(proteins$gene_name, liana_genes)
log_info(glue("Number of genes to keep: {length(genes_to_keep)}"))
genes_input <- proteins %>%
    select(gene_name, uniprot, hgnc_symbol, ensembl) %>%
    filter(gene_name %in% genes_to_keep) %>%
    distinct()

interactions_to_keep <- rowSums(apply(liana_db %>%
    select(starts_with("ligand_subunit") | starts_with("receptor_subunit")), 2, function(col) {
    col %in% c(missing_genes)
})) == 0
liana_db <- liana_db[interactions_to_keep, ]

log_info(glue("Number of interactions: {nrow(liana_db)}"))

# Protein input
# mandatory: uniprot, protein_name
protein_input <- proteins %>%
    select(uniprot, protein_name) %>%
    distinct()

# Complex input
# mandatory: complex_name, complex_subunits
ligand_complexes <- liana_db %>%
    filter(ligand_subunit_2 != "" | is.na(ligand_subunit_2)) %>%
    separate(source_protein, paste0("uniprot_", seq_len(5)), sep = "_", remove = FALSE) %>%
    select(source_genesymbol, starts_with("uniprot_")) %>%
    rename(complex_name = source_genesymbol)
receptor_complexes <- liana_db %>%
    filter(receptor_subunit_2 != "" | is.na(receptor_subunit_2)) %>%
    separate(target_protein, paste0("uniprot_", seq_len(5)), sep = "_", remove = FALSE) %>%
    select(target_genesymbol, starts_with("uniprot_")) %>%
    rename(complex_name = target_genesymbol)

complex_input <- rbind(ligand_complexes, receptor_complexes) %>% distinct()
complex_input[is.na(complex_input)] <- ""
complex_uniprot <- complex_input %>% filter(!is.na(uniprot_2))
log_info(glue("Number of complexes: {nrow(complex_uniprot)}"))

# Interaction Input
# mandatory:  “partner_a”; “partner_b”; “annotation_strategy”; “source”
interactions_input <- liana_db %>%
    select(source_protein, target_protein) %>%
    rename(partner_a = source_protein, partner_b = target_protein) %>%
    mutate(annotation_strategy = "user_curated", version = "CellPhoneDBcore4.1", sources = "User curated")


liana_proteins <- liana_db %>%
    select(source_protein, target_protein) %>%
    unlist() %>%
    unique() %>%
    str_split(., "_") %>%
    unlist() %>%
    unique()
log_info(glue("Number of proteins: {length(liana_proteins)}"))
missing_proteins <- setdiff(liana_proteins, protein_input$uniprot)
log_info(glue("Number of missing proteins: {length(missing_proteins)}"))

# missing cols
missing_cols <- setdiff(gene_input_cols, colnames(genes_input))
genes_input[missing_cols] <- ""
missing_cols <- setdiff(protein_input_cols, colnames(protein_input))
protein_input[missing_cols] <- ""
missing_cols <- setdiff(complex_input_cols, colnames(complex_input))
complex_input[missing_cols] <- ""
complex_input$version <- "CellPhoneDBcore4.1"
missing_cols <- setdiff(interaction_input_cols, colnames(interactions_input))
interactions_input[missing_cols] <- ""
log_info(glue("Number of interactions: {nrow(interactions_input)}"))

log_info("Save updated files...")
write.csv(genes_input, glue("{args$output_dir}/gene_input.csv"), row.names = FALSE)
write.csv(protein_input, glue("{args$output_dir}/protein_input.csv"), row.names = FALSE)
write.csv(complex_input, glue("{args$output_dir}/complex_input.csv"), row.names = FALSE)
write.csv(interactions_input, glue("{args$output_dir}/interaction_input.csv"), row.names = FALSE)

log_info("COMPLETED!")
