#' Setup default argparser
#'
#' Set up a default argument parser, with two default arguments: log_level and output_dir
#'
#' @param description Description for script
#' @param output_default Default output directory
#' @return parser object
#'
#' @examples setup_default_argparser(description = "Example")
#' @export
setup_default_argparser <- function(description = "", default_output = "output") {
  parser <- argparse::ArgumentParser(
    description = description, python_cmd = NULL
  )
  parser$add_argument("-ll", "--log_level",
    type = "integer",
    default = "4",
    help = "Log level: 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, 5=DEBUG"
  )
  parser$add_argument("-o", "--output_dir",
    type = "character",
    default = default_output, help = "Directory to save output"
  )
  return(parser)
}
