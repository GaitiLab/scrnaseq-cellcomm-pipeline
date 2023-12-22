#' Get file information
#'
#' Split the filename based on double underscores to extract several characteristics of the file.
#'
#' @param p path to file
#' @return file_info in array
#' @export
#' @examples get_file_info("cell2cell__6245__Tumour_edge__0.csv")
#' @importFrom stringr str_split
get_file_info <- function(p) {
    file_info <- str_split(get_name(p), "__", simplify = TRUE)
    return(file_info)
}

#' Get metadata for Cell2Cell
#' @param paths paths to files
#' @param files_info_colnames column names for metadata
#' @return metadata dataframe
#' @export
get_metadata <- function(paths, file_info_colnames = c("method", "sample_id", "region", "iteration")) {
    # Metadata for sample
    files_info <- data.frame(do.call(
        rbind,
        lapply(paths, get_file_info)
    ))
    colnames(files_info) <- file_info_colnames
    return(files_info)
}

#' Get sample IDs
#'
#' @param metadata_df metadata dataframe
#' @return sample_ids in array
#' @export
#'
#' @importFrom dplyr select pull %>%
get_sample_ids <- function(metadata_df) {
    # Get unique sample IDs
    return(metadata_df %>%
        select(sample_id) %>%
        unique() %>%
        pull())
}
#' Format output of CellChat
#'
#' @param p path to file
#' @return obj formatted output
#' @export
#' @importFrom stringr str_split
#' @importFrom dplyr mutate %>% select
#' @importFrom glue glue
format_cellchat <- function(p, file_info_colnames = c("method", "sample_id", "region", "iteration")) {
    file_info <- c(get_file_info(p))
    names(file_info) <- file_info_colnames
    obj <- readRDS(p) %>%
        mutate(
            method = file_info["method"],
            sample_id = file_info["sample_id"],
            region = file_info["region"],
            iteration = glue("pval_{file_info['iteration']}"),
            pval = as.numeric(pval)
        )
    return(obj)
}
#' Format output of Liana
#'
#' @param p path to file
#' @return obj formatted output
#' @export
#'
#' @importFrom liana liana_aggregate
#' @importFrom dplyr mutate rename %>% select rename
#' @importFrom glue glue
format_liana <- function(p, file_info_colnames = c("method", "sample_id", "region", "iteration")) {
    file_info <- get_file_info(p)
    names(file_info) <- file_info_colnames
    obj <- readRDS(p) %>%
        liana_aggregate() %>%
        dplyr::rename(pval = aggregate_rank) %>%
        select(source, target, ligand.complex, receptor.complex, pval) %>%
        mutate(
            method = file_info["method"],
            sample_id = file_info["sample_id"],
            region = file_info["region"],
            iteration = glue("pval_{file_info['iteration']}"),
            pval = as.numeric(pval),
            ligand.complex = str_replace_all(ligand.complex, "_", ":"), receptor.complex = str_replace_all(receptor.complex, "_", ":"), interaction = paste0(ligand.complex, "__", receptor.complex),
        )
    return(obj)
}

#' Format output of Cell2Cell
#'
#' @param p path to file
#' @return obj format
#' @export
#'
#' @importFrom data.table fread
#' @importFrom glue glue
#' @importFrom reshape2 melt
#' @importFrom stringr str_split str_replace_all
#' @importFrom dplyr mutate %>% rename
#' @importFrom data.table fread
#'
format_cell2cell <- function(p, file_info_colnames = c("method", "sample_id", "region", "iteration")) {
    #  ---- Constants ---- #
    str_format_source_target <- c(";" = "__")
    str_format_interaction <- c(
        "\\(" = "", "\\)" = "", # remove brackets
        "_" = ":", # subunits
        ", " = "__", # L-R
        "\'" = "" # remove quotes (redundant)
    )

    file_info <- c(get_file_info(p))
    names(file_info) <- file_info_colnames
    obj <- data.frame(fread(p, header = TRUE, sep = "\t"),
        check.names = FALSE
    ) %>%
        rename(interaction = V1) %>%
        reshape2::melt(id.vars = c("interaction"), variable.name = "source_target", value.name = "pval") %>%
        mutate(
            method = file_info["method"],
            sample_id = file_info["sample_id"],
            region = file_info["region"],
            iteration = glue("pval_{file_info['iteration']}"),
            pval = as.numeric(pval)
        ) %>%
        rowwise() %>%
        mutate(
            interaction = str_replace_all(
                interaction,
                str_format_interaction
            ),
            source_target = str_replace_all(
                as.character(source_target),
                str_format_source_target
            )
        )
    return(obj)
}
