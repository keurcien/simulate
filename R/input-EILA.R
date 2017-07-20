#' Input tools
#'
#' \code{create_input_EILA} creates an object storing the input matrices 
#' required by \code{EILA::eila}.
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' @param Hadm a haplotype matrix.
#' @param position a numerical vector.
#' 
#' @return 
#'
#' @export
#'
create_input_EILA = function(H1, 
                             H2, 
                             Hadm,
                             position){
  obj.eila <- NULL
  obj.eila$anc1 <- haplo_to_geno(H1)  
  obj.eila$anc2 <- haplo_to_geno(H2)
  obj.eila$admixed <- haplo_to_geno(Hadm)
  obj.eila$position <- position
  return(obj.eila)
}

#' Input tools
#'
#' \code{EILA_from_pcadapt} creates an object storing the input matrices 
#' required by \code{EILA::eila}.
#'
#' @param geno a genotype matrix. 
#' @param pop a vector of integers.
#' @param ancestral.1 an integer.
#' @param ancestral.2 an integer.
#' @param admixed an integer.
#' @param position a numerical vector.
#' 
#' @return 
#'
#' @export
#'
EILA_from_pcadapt = function(geno, 
                             pop, 
                             ancestral.1, 
                             ancestral.2, 
                             admixed, 
                             position){
  obj.eila <- NULL
  obj.eila$anc1 <- geno[, pop == ancestral.1]
  obj.eila$anc2 <- geno[, pop == ancestral.2]
  obj.eila$admixed <- geno[, pop == admixed]
  obj.eila$position <- position
  return(obj.eila)
}