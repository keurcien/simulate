#' Sample introgression regions
#'
#' \code{sample_introgression_regions} returns a list of introgressed alleles.
#'
#' @param nSNP an integer.
#' @param n.windows an integer.
#' @param n.regions an integer.
#' 
#' @return 
#'
#' @export
#'
sample_introgression_regions = function(nSNP, n.windows, n.regions){
  idx <- sort(unique(sample(1:n.windows, size = n.regions)))
  window.size <- nSNP / n.windows
  gt.df <- data.frame(start = (idx - 1) * window.size + 1 , 
                      end = idx * window.size)
  introgression.regions <- NULL
  for (k in 1:n.regions){
    introgression.regions <- c(introgression.regions, 
                               (gt.df$start[k]):(gt.df$end[k]))
  }
  return(list(gt.df = gt.df, 
              introgression.regions = introgression.regions))
}