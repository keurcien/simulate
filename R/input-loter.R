#' Input tools
#'
#' \code{create_input_loter} 
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' @param Hadm a haplotype matrix.
#' @param loter.output a character string.
#' 
#' @return 
#'
#' @export
#'
create_input_loter = function(H1, 
                              H2, 
                              Hadm, 
                              loter.output){
  
  G1 <- haplo_to_geno(H1)  
  G2 <- haplo_to_geno(H2)
  Gadm <- haplo_to_geno(Hadm)
  
  write.table(t(G1), 
              paste0(loter.output, "loter_refpop1_geno.txt"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  write.table(t(G2), 
              paste0(loter.output, "loter_refpop2_geno.txt"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  write.table(t(Gadm), 
              paste0(loter.output, "loter_admpop_geno.txt"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  write.table(t(H1), 
              paste0(loter.output, "loter_refpop1_haplo.txt"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  write.table(t(H2), 
              paste0(loter.output, "loter_refpop2_haplo.txt"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  write.table(t(Hadm), 
              paste0(loter.output, "loter_admpop_haplo.txt"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
}