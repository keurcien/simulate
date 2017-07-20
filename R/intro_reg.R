#' Sample introgression regions
#'
#' \code{sample_intro_reg} returns a list of introgressed alleles.
#'
#' @param nSNP an integer.
#' @param n.reg an integer.
#' @param intro.size an integer.
#' 
#' @return 
#'
#' @export
#'
sample_intro_reg = function(nSNP, n.reg = 5, intro.size){
  ### Introgression region
  intro.reg <- NULL
  idx <- sample(1:nSNP, size = n.reg)
  gt.df <- data.frame(beg = pmax(1, idx - intro.size), 
                      end = pmin(nSNP, idx + intro.size))
  for (k in 1:n.reg){
    intro.reg <- c(intro.reg, (gt.df$beg[k]):(gt.df$end[k]))
  }
  return(list(gt.df = gt.df, intro.reg = intro.reg))
}

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
  idx <- sample(1:n.windows, size = n.regions)
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