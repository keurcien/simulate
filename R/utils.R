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

#' Utils
#'
#' \code{stat.by.window} split a vector of statistics into a certain number
#' windows.
#'
#' @param stat a numeric vector.
#' @param n.windows an integer specifying the number of windows.

#' 
#' @return 
#'
#' @export
#'
stat.by.window = function(stat, 
                          n.windows, 
                          window.size){
  sbw.max <- vector(mode = "numeric", n.windows)
  sbw.mean <- vector(mode = "numeric", n.windows)
  subset.min <- (1:n.windows - 1) * window.size + 1
  subset.max <- (1:n.windows) * window.size
  sbw.max <- sapply(1:n.windows, FUN = function(h){max(stat[subset.min[h]:subset.max[h]], na.rm = TRUE)})
  sbw.mean <- sapply(1:n.windows, FUN = function(h){mean(stat[subset.min[h]:subset.max[h]], na.rm = TRUE)})
  return(list(stat.max = sbw.max, stat.mean = sbw.mean))
}

#' Utils
#'
#' \code{top} 
#'
#' @param idx a numeric vector.
#' @param gt a numeric vector.
#' @param n an integer.
#' 
#' @return 
#'
#' @export
#'
top = function(idx, gt, n){
  sum(idx[1:n] %in% gt)
}