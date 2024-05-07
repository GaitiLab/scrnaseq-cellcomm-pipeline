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
        description = "Update CellChat database w/ LIANA",
    )
    parser$add_argument("--liana_db",
        default = "",
        type = "character", help = "Path to LIANA database"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/001_data_local/interactions_db_v2/")
    args$liana_db <- "001_data_local/interactions_db_v2/liana_db.rds"
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
pacman::p_load_gh("jinworks/CellChat")

# Ref: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Update-CellChatDB.html

log_info("Load CellChat database...")
cellchat_db <- CellChatDB.human
cellchat_db_interaction <- cellchat_db$interaction
cellchat_db_complex <- cellchat_db$complex
cellchat_db_cofactor <- cellchat_db$cofactor
cellchat_db_geneinfo <- cellchat_db$geneInfo

log_info("Load LIANA database...")
liana_db <- readRDS(ifelse(file.exists(args$liana_db),
    args$liana_db, glue("{args$output_dir}/liana_db.rds")
))
log_info(glue("Number of interactions: {nrow(liana_db)}"))

# Create new CellChatDB
log_info("Create new CellChatDB...")
cellchatDB_Omni <- liana:::cellchat_formatDB(
    ccDB = CellChat::CellChatDB.human,
    op_resource = liana_db,
    exclude_anns = c()
)
log_info("Save CellChatDB...")
saveRDS(cellchatDB_Omni, file = glue("{args$output_dir}/cellchat_db.rds"))
log_info("COMPLETED!")
