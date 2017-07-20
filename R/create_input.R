#' Input tools
#'
#' \code{create_all_inputs} 
#'
#' @param refpop1_integer a haplotype matrix. 
#' @param refpop2_integer a haplotype matrix.
#' @param refpop1 a haplotype matrix. 
#' @param refpop2 a haplotype matrix. 
#' @param H.integer a haplotype matrix.
#' @param H a haplotype matrix.
#' @param true.ancestry.matrix a matrix.
#' @param directory a character string.
#' @param pos_map a numeric vector.
#' @param gen_map a numeric vector.
#' @param chr an integer.
#' @param ref a numeric vector.
#' @param alt a numeric vector.
#' @param pop a numeric vector.
#' @param parameters a data.frame.
#' @param gt.df a data.frame.
#' 
#' @return 
#'
#' @export
#'
create_all_inputs = function(refpop1_integer, 
                             refpop2_integer,
                             refpop1,
                             refpop2,
                             H.integer,
                             H,
                             true.ancestry.matrix,
                             directory, 
                             pos_map, 
                             gen_map, 
                             chr,
                             ref,
                             alt,
                             pop, 
                             parameters, 
                             gt.df){
  
  ### pcadapt ###
  input.pcadapt <- create_input_pcadapt(H1 = refpop1_integer, 
                                        H2 = refpop2_integer,
                                        Hadm = H.integer)
  write.table(input.pcadapt, 
              paste0(directory, "/simu.pcadapt"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  
  ### RFMix ###
  create_input_RFMix(H1 = refpop1_integer, 
                     H2 = refpop2_integer,
                     Hadm = H.integer, 
                     genetic.map = gen_map, 
                     RFMix.output = paste0(directory, "/RFMix"))
  
  ### loter ###
  create_input_loter(H1 = H1.integer, 
                     H2 = H2.integer, 
                     Hadm = H.integer, 
                     paste0(dir.name, "/"))  
  
  ### Hapmix ###
  create_input_Hapmix(H1 = H1.integer, 
                      H2 = H2.integer, 
                      Gadm = haplo_to_geno(H.integer),
                      gen_map = gen_map,
                      phys_map = pos_map,
                      chr = chr, 
                      Hapmix.output = directory)
  
  write.table(pop, 
              paste0(directory, "/pop.txt"), 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
  
  write.table(parameters, 
              paste0(directory, "/parameters.txt"), 
              col.names = TRUE, 
              row.names = FALSE,
              quote = FALSE)
  
  write.table(gt.df$gt.df, 
              paste0(directory, "/gt.txt"), 
              col.names = TRUE, 
              row.names = FALSE, 
              quote = FALSE)
  
  write.table(true.ancestry.matrix, 
              paste0(directory, "/adm_true_ancestry.txt"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  
  write.table(H, 
              paste0(directory, "/admpop.phased"),
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
}