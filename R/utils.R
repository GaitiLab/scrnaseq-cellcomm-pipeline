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