#' Input tools
#'
#' \code{create_input_Hapmix} generates input files for Hapmix.
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' @param Gadm a genotype matrix.
#' @param gen_map a numerical vector.
#' @param phys_map a numerical vector.
#' @param chr an integer.
#' @param Hapmix.output a character string.
#' @param ancestral.1 a character string.
#' @param ancestral.2 a character string.
#' @param admixed a character string.
#' @param rates.digits an integer.
#' 
#' @importFrom utils write.table
#' @importFrom LEA lfmm2geno
#' 
#' @return 
#'
#' @export
#'
create_input_Hapmix = function(H1, 
                               H2, 
                               Gadm, 
                               gen_map, 
                               phys_map, 
                               chr, 
                               Hapmix.output,
                               ancestral.1 = "ANC1", 
                               ancestral.2 = "ANC2", 
                               admixed = "AA", 
                               rates.digits = 8){
  cat("Creating the genotype file for ancestral population 1...")
  write.table(t(H1), 
              paste0(Hapmix.output, "/", ancestral.1, ".lfmm"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(Hapmix.output, "/", ancestral.1, ".lfmm"), 
                 paste0(Hapmix.output, "/", ancestral.1, ".geno"), 
                 force = TRUE)
  file.remove(paste0(Hapmix.output, "/", ancestral.1, ".lfmm")) 
  cat("Creating the genotype file for ancestral population 2...")
  write.table(t(H2), 
              paste0(Hapmix.output, "/", ancestral.2, ".lfmm"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(Hapmix.output, "/", ancestral.2, ".lfmm"), 
                 paste0(Hapmix.output, "/", ancestral.2, ".geno"), 
                 force = TRUE)
  file.remove(paste0(Hapmix.output, "/", ancestral.2, ".lfmm"))
  cat("Creating the genotype file for the admixed population...")
  write.table(t(Gadm), 
              paste0(Hapmix.output, "/", admixed, ".lfmm"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(Hapmix.output, "/", admixed, ".lfmm"), 
                 paste0(Hapmix.output, "/", admixed, ".geno"), 
                 force = TRUE)
  file.remove(paste0(Hapmix.output, "/", admixed, ".lfmm"))
  dt <- array(dim = c(nrow(H1), ncol = 4))
  dt[, 1] <- paste0("rs", 1:nrow(H1))
  dt[, 2] <- chr
  dt[, 3] <- format(round(gen_map * 1e-2, rates.digits), nsmall = rates.digits)
  dt[, 4] <- phys_map
  cat("Creating the snpfiles...")
  write.table(dt, 
              paste0(Hapmix.output, "/", ancestral.1, "snpfile"), 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
  write.table(dt, 
              paste0(Hapmix.output, "/", ancestral.2, "snpfile"), 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
  write.table(dt, 
              paste0(Hapmix.output, "/", admixed, "snpfile"), 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
  dt.ind <- as.data.frame(matrix(nrow = ncol(Gadm), ncol = 3))
  dt.ind[[1]] <- paste0(rep("IND_", ncol(Gadm)), 1:ncol(Gadm))
  dt.ind[[2]] <- "F"
  dt.ind[[3]] <- "Hybrid"
  write.table(dt.ind, 
              paste0(Hapmix.output, "/", admixed, ".ind"), 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
  cat("DONE\n")
  cat("Creating the rates file...")
  l1 <- paste0(":sites:", nrow(Gadm))
  l2 <- stringr::str_c(phys_map, collapse = " ")
  l3 <- stringr::str_c(format(x = gen_map * 1e-2, 
                              scientific = FALSE, 
                              nsmall = rates.digits), 
                       collapse = " ")
  con <- file(paste0(Hapmix.output, "/rates"))
  writeLines(paste0(l1, "\n", l2, "\n", l3), con = con)
  close(con = con)
  cat("DONE\n")
}
