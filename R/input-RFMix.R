#' RFMix
#'
#' \code{create_input_RFMix}
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' @param Hadm a haplotype matrix.
#' @param genetic.map a numerical vector.
#' @param RFMix.output a character string.
#' 
#' @return 
#'
#' @export
#'
create_input_RFMix = function(H1, 
                              H2, 
                              Hadm, 
                              genetic.map, 
                              RFMix.output = "RFMix"){
  G <- cbind(H1, H2, Hadm)
  write.table(t(G), 
              paste0(RFMix.output, "_alleles.lfmm"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  LEA::lfmm2geno(paste0(RFMix.output, "_alleles.lfmm"), 
                 paste0(RFMix.output, "_alleles.geno"), 
                 force = TRUE)
  file.remove(paste0(RFMix.output, "_alleles.lfmm"))
  file.rename(paste0(RFMix.output, "_alleles.geno"), 
              paste0(RFMix.output, "_alleles.txt"))
  formatted.genetic.map <- format(genetic.map, 
                                  scientific = FALSE, 
                                  nsmall = 8) 
  write.table(formatted.genetic.map, 
              paste0(RFMix.output, "_markerLocation.txt"), 
              sep = "\n", 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
  
  pop <- c(rep(1, ncol(H1)), rep(2, ncol(H2)), rep(0, ncol(Hadm)))
  l3 <- stringr::str_c(pop, collapse = " ")
  con <- file(paste0(RFMix.output, "_classes.txt"))
  writeLines(l3, con = con)
  close(con = con)
}