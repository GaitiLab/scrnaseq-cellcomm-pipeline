#' Extract name from filepath
#'
#' @param filepath filepath
#' @return name of file without extension
#'
#' @examples examples get_name("R/my_script.R") returns 'my_script'
#' @export
#' @importFrom tools file_path_sans_ext
#' @importFrom fs path_file
#'
get_name <- function(filepath) {
    return(file_path_sans_ext(path_file(filepath)))
}

#' Create directory
#'
#' @description Create directory if it does not exist
#'
#' @param dir_path Path to directory to be created
#'
#' @examples create_dir("/Users/johndoe/my_new_dir")
#' @export
#' @importFrom glue glue
create_dir <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        log_info(glue("Creating directory {dir_path}"))
        dir.create(dir_path, recursive = TRUE)
    } else {
        log_warn(glue("Directory {dir_path} already exists."))
    }
}

#' Extract patient id from filepath
#'
#' @description Extract patient id from filepath from following format:
#' liana__6237_2222190_A__Cortex__Batch_1, where 6237_2222190 is patient id
#'
#' @param filepath filepath
#' @return patient id
#'
#' @export
#' @importFrom stringr str_split
extract_patient_id <- function(filepath, platform = "parsebio") {
    if (platform == "parsebio") {
    return(paste(stringr::str_split(get_name(filepath),
        "_",
        simplify = TRUE
    )[c(3, 4)], collapse = "_"))
    } else if (platform == "10x") {
        return(stringr::str_split(stringr::str_split(get_name(files[1]), "__", simplify = TRUE)[2], "_", simplify = TRUE)[1])
    }
}

#' Obtain current date
#' @description Obtain current date in format YYYYMMDD
#' @return current date
#' @export
get_current_date <- function() {
    return(format(Sys.time(), "%Y%m%d"))
}


#' Combine p-values
#' @description Combine p-values using Fisher's method
#' @param p_values vector of p-values
#' @return combined p-value
#' @export
#' @importFrom metap sumlog
get_p <- function(x) {
    if (length(x) > 1) {
        return(sumlog(x)$p)
    } else {
        return(x)
    }
}