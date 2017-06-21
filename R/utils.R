#' Utils
#'
#' \code{number.to.string} converts a numeric value to a character string.
#'
#' @param n a numeric value.
#' 
#' @return 
#'
#' @export
#'
number.to.string = function(n){
  str <- format(n, scientific = FALSE)
  tokens <- strsplit(str, split = "")
  figures.only <- tokens[[1]][tokens[[1]] != "."]
  return(paste(figures.only, collapse = ""))
}