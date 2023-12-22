#' Set up logging
#'
#' @param log_level level of logging (1-5)
#' @param log_file path to log file (optional)
#' @return logger object
#'
#' @examples logr <- init_logging(3)
#' @export
init_logging <- function(log_level, log_file = NULL) {
  log_level_options <- c(
    `1` = "FATAL", `2` = "ERROR", `3` = "WARN", `4` = "INFO",
    `5` = "DEBUG"
  )
  if (!is.null(log_file)) {
    console_appender <- log4r::console_appender(layout = log4r::default_log_layout())
    file_appender <- log4r::file_appender(log_file,
      append = FALSE,
      layout = log4r::default_log_layout()
    )
    return(log4r::logger(
      threshold = log_level_options[as.character(log_level)],
      appenders = list(console_appender, file_appender)
    ))
  }
  return(log4r::logger(
    threshold = log_level_options[as.character(log_level)],
    appenders = log4r::console_appender(layout = log4r::default_log_layout())
  ))
}

#' Logging functions: log_info
#'
#' @param logr logger object
#' @param ... message
#'
#' @examples log_info("Hello world!")
#' @export
log_info <- function(...) {
  log4r::info(logr, paste0(...))
}

#' Logging functions: log_error
#'
#' @param logr logger object
#' @param ... message
#'
#' @examples log_error("Hello world!")
#'
#' @export
log_error <- function(...) {
  log4r::error(logr, paste0(...))
}

#' Logging functions: log_fatal
#'
#' @param logr logger object
#' @param ... message
#'
#' @examples log_fatal("Hello world!")
#'
#' @export
log_fatal <- function(...) {
  log4r::fatal(logr, paste0(...))
}

#' Logging functions: log_debug
#'
#' @param logr logger object
#' @param ... message

#' @examples log_debug("Hello world!")
#'
#' @export
log_debug <- function(...) {
  log4r::debug(logr, paste0(...))
}

#' Logging functions: log.warn
#'
#' @param logr logging object
#' @param ... message
#'
#' @examples log_warn("Hello world!")
#' @export
log_warn <- function(...) {
  log4r::warn(logr, paste0(...))
}
