#' Input tools
#'
#' \code{create_input_pcadapt} creates a genotype matrix in the \code{pcadapt} 
#' format.
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' @param Hadm a haplotype matrix.
#' 
#' @return 
#'
#' @export
#'
create_input_pcadapt = function(H1, H2, Hadm){
  nIND <- (ncol(H1) + ncol(H2) + ncol(Hadm)) / 2
  Gmat <- matrix(0, nrow = nrow(Hadm), ncol = nIND)
  Gmat[, 1:(ncol(H1) / 2)] <- haplo_to_geno(H1)
  Gmat[, (ncol(H1) / 2 + 1):(ncol(H1) / 2 + ncol(H2) / 2)] <- haplo_to_geno(H2)
  Gmat[, (ncol(H1) / 2 + ncol(H2) / 2 + 1):nIND] <- haplo_to_geno(Hadm)
  return(Gmat)
}