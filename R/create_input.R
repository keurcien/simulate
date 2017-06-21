#' Input tools
#'
#' \code{create_input_pcadapt} creates a genotype matrix in the \code{pcadapt} 
#' format.
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' @param H a haplotype matrix.
#' 
#' @return 
#'
#' @export
#'
create_input_pcadapt = function(H1, H2, H){
  G1 <- haplo_to_geno(H1)  
  G2 <- haplo_to_geno(H2)
  G <- haplo_to_geno(H)
  nIND <- ncol(G1) + ncol(G2) + ncol(G)
  Gmat <- matrix(0, nrow = nrow(G), ncol = nIND)
  Gmat[, 1:ncol(G1)] <- G1
  Gmat[, (ncol(G1) + 1):(ncol(G1) + ncol(G2))] <- G2
  Gmat[, (ncol(G1) + ncol(G2) + 1):(ncol(G1) + ncol(G2) + ncol(G))] <- G
  return(Gmat)
}

#' Input tools
#'
#' \code{create_input_hapmix} generates input files for Hapmix.
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' @param H a haplotype matrix.
#' @param gen_map a numerical vector.
#' @param phys_map a numerical vector.
#' @param chr an integer.
#' @param ancstrl.1 a character string.
#' @param ancstrl.2 a character string.
#' @param admxd a character string.
#' @param rates.digits an integer.
#' 
#' @importFrom utils write.table
#' @importFrom LEA lfmm2geno
#' 
#' @return 
#'
#' @export
#'
create_input_hapmix = function(H1, H2, H, gen_map, phys_map, chr, 
                               ancstrl.1 = "ANC1", ancstrl.2 = "ANC2", 
                               admxd = "AA", rates.digits = 8){
  cat("Creating the genotype file for ancestral population 1...")
  write.table(t(H1), paste0(ancstrl.1, ".lfmm"), 
              col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(ancstrl.1, ".lfmm"), paste0(ancstrl.1, ".geno"), 
                 force = TRUE)
  cat("Creating the genotype file for ancestral population 2...")
  write.table(t(H2), paste0(ancstrl.2, ".lfmm"), 
              col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(ancstrl.2, ".lfmm"), paste0(ancstrl.2, ".geno"), 
                 force = TRUE)
  cat("Creating the genotype file for the admixed population...")
  G <- haplo_to_geno(H)
  write.table(t(G), paste0(admxd, ".lfmm"), 
              col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(admxd, ".lfmm"), paste0(admxd, ".geno"), 
                 force = TRUE)
  
  dt <- array(dim = c(nrow(H1), ncol = 4))
  dt[, 1] <- paste0("rs", 1:nrow(H1))
  dt[, 2] <- chr
  dt[, 3] <- format(round(gen_map, rates.digits), nsmall = rates.digits)
  dt[, 4] <- phys_map
  cat("Creating the snpfiles...")
  write.table(dt, paste0(ancstrl.1, "snpfile"), col.names = FALSE, 
              row.names = FALSE, quote = FALSE)
  write.table(dt, paste0(ancstrl.2, "snpfile"), col.names = FALSE, 
              row.names = FALSE, quote = FALSE)
  write.table(dt, paste0(admxd, "snpfile"), col.names = FALSE, 
              row.names = FALSE, quote = FALSE)
  dt.ind <- as.data.frame(matrix(nrow = ncol(G), ncol = 3))
  dt.ind[[1]] <- paste0(rep("IND_", ncol(G)), 1:ncol(G))
  dt.ind[[2]] <- "F"
  dt.ind[[3]] <- "Hybrid"
  write.table(dt.ind, paste0(admxd, ".ind"), col.names = FALSE, 
              row.names = FALSE, quote = FALSE)
  cat("DONE\n")
  cat("Creating the rates file...")
  l1 <- paste0(":sites:", nrow(G))
  l2 <- stringr::str_c(phys_map, collapse = " ")
  l3 <- stringr::str_c(format(gen_map, scientific = FALSE, 
                              nsmall = rates.digits), 
                       collapse = " ")
  con <- file("rates")
  writeLines(paste0(l1, "\n", l2, "\n", l3), con = con)
  close(con = con)
  cat("DONE\n")
}

#' Input tools
#'
#' \code{create_input_eila} creates an object storing the input matrices 
#' required by \code{EILA::eila}.
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' @param H a haplotype matrix.
#' @param position a numerical vector.
#' 
#' @return 
#'
#' @export
#'
create_input_eila = function(H1, H2, H, position){
  obj.eila <- NULL
  obj.eila$anc1 <- haplo_to_geno(H1)  
  obj.eila$anc2 <- haplo_to_geno(H2)
  obj.eila$admixed <- haplo_to_geno(H)
  obj.eila$position <- position
  return(obj.eila)
}

#' Input tools
#'
#' \code{eila_from_pcadapt} creates an object storing the input matrices 
#' required by \code{EILA::eila}.
#'
#' @param geno a genotype matrix. 
#' @param pop a vector of integers.
#' @param anc1 an integer.
#' @param anc2 an integer.
#' @param admixed an integer.
#' @param position a numerical vector.
#' 
#' @return 
#'
#' @export
#'
eila_from_pcadapt = function(geno, pop, anc1, anc2, admixed, position){
  obj.eila <- NULL
  obj.eila$anc1 <- geno[, pop == anc1]
  obj.eila$anc2 <- geno[, pop == anc2]
  obj.eila$admixed <- geno[, pop == admixed]
  obj.eila$position <- position
  return(obj.eila)
}

#' RFMix
#'
#' \code{create_input_rfmix}
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' @param H3 a haplotype matrix.
#' @param gen_map a numerical vector.
#' @param rfmix.output a character string.
#' 
#' @return 
#'
#' @export
#'
create_input_rfmix = function(H1, H2, H3, gen_map, rfmix.output = "rfmix"){
  G <- cbind(H1, H2, H3)
  write.table(t(G), paste0(rfmix.output, "_alleles.lfmm"), 
              col.names = FALSE, row.names = FALSE)
  LEA::lfmm2geno(paste0(rfmix.output, "_alleles.lfmm"), 
                 paste0(rfmix.output, "_alleles.geno"), force = TRUE)
  file.remove(paste0(rfmix.output, "_alleles.lfmm"))
  file.rename(paste0(rfmix.output, "_alleles.geno"), 
              paste0(rfmix.output, "_alleles.txt"))
  formatted.gen_map <- format(gen_map, scientific = FALSE, nsmall = 8) 
  write.table(formatted.gen_map, paste0(rfmix.output, "_markerLocation.txt"), 
              sep = "\n", col.names = FALSE, row.names = FALSE, quote = FALSE)
  pop <- c(rep(1, ncol(H1)), rep(2, ncol(H2)), rep(0, ncol(H3)))
  l3 <- stringr::str_c(pop, collapse = " ")
  con <- file(paste0(rfmix.output, "_classes.txt"))
  writeLines(l3, con = con)
  close(con = con)
}

create_input_loter = function(H1, H2, H, loter.output){
  G1 <- haplo_to_geno(H1)  
  G2 <- haplo_to_geno(H2)
  G <- haplo_to_geno(H)
  write.table(t(G1), paste0(loter.output, "G1_loter.txt"), col.names = FALSE, row.names = FALSE)
  write.table(t(G2), paste0(loter.output, "G2_loter.txt"), col.names = FALSE, row.names = FALSE)
  write.table(t(G), paste0(loter.output, "G_loter.txt"), col.names = FALSE, row.names = FALSE)
  write.table(t(H1), paste0(loter.output, "H1_loter.txt"), col.names = FALSE, row.names = FALSE)
  write.table(t(H2), paste0(loter.output, "H2_loter.txt"), col.names = FALSE, row.names = FALSE)
  write.table(t(H), paste0(loter.output, "H_loter.txt"), col.names = FALSE, row.names = FALSE)
}

create_all_inputs = function(H1, H2, H, directory, pos_map, gen_map, pop, df, gt.df){
  input.pcadapt <- create_input_pcadapt(H1, H2, H)
  write.table(input.pcadapt, paste0(directory, "/simu.pcadapt"), col.names = FALSE, row.names = FALSE)
  input.eila <- create_input_eila(H1, H2, H, position = pos_map)
  create_input_rfmix(H1, H2, H, gen_map = gen_map, rfmix.output = paste0(directory, "/rfmix"))
  #create_input_loter(H1, H2, H3, paste0(dir.name, "/"))  
  write.table(pop, paste0(directory, "/pop.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(df, paste0(directory, "/parameters.txt"), col.names = TRUE, quote = FALSE)
  write.table(gt.df, paste0(directory, "/gt.txt"), col.names = TRUE, quote = FALSE)
}