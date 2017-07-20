#' Input tools
#'
#' \code{create_input_RASPberry} generates input files for RASPberry.
#'
#' @param refpop1_phased a haplotype matrix. 
#' @param refpop2_phased a haplotype matrix.
#' @param H a genotype matrix.
#' @param pos_map a numerical vector.
#' @param ref a string vector.
#' @param alt a string vector.
#' @param chr an integer.
#' 
#' @export
#'
create_input_RASPberry = function(refpop1_phased,
                                  refpop2_phased, 
                                  H,
                                  pos_map,
                                  ref,
                                  alt,
                                  chr){
  nSNP <- nrow(refpop1_phased)
  nHAP <- ncol(refpop2_phased)
  subseq_A <- paste0("IND", 1:(nHAP / 2), "_A")
  subseq_B <- paste0("IND", 1:(nHAP / 2), "_B")
  subseq <- vector(mode = "character", length = nHAP)
  subseq[seq(1, nHAP, by = 2)] <- subseq_A
  subseq[seq(2, nHAP, by = 2)] <- subseq_B
  
  H1.table <- cbind(data.frame(rsID = paste0("rs", 1:nSNP + 1e7), 
                               position = pos_map), 
                    refpop1_phased)
  H2.table <- cbind(data.frame(rsID = paste0("rs", 1:nSNP + 1e7), 
                               position = pos_map), 
                    refpop2_phased)
  colnames(H1.table) <- c("rsID", "position", subseq)
  colnames(H2.table) <- c("rsID", "position", subseq)
  bim.table <- data.frame(Chr = rep(chr, nSNP), 
                          rsID = H1.table$rsID,
                          unknown = rep(0, nSNP),
                          position = H1.table$position,
                          REF = ref,
                          ALT = alt)
  map.table <- data.frame(Chr = rep(chr, nSNP), 
                          rsID = H1.table$rsID,
                          unknown = rep(0, nSNP),
                          position = H1.table$position)
  Hadm <- matrix(0, nrow = 2 * nSNP, ncol = (nHAP / 2))
  sbs.snp.1 <- seq(1, 2 * nSNP, by = 2)
  sbs.snp.2 <- seq(2, 2 * nSNP, by = 2)
  sbs.hap.1 <- seq(1, nHAP, by = 2)
  sbs.hap.2 <- seq(2, nHAP, by = 2)
  Hadm[sbs.snp.1, ] <-H[, sbs.hap.1]
  Hadm[sbs.snp.2, ] <-H[, sbs.hap.2]
  nIND <- nHAP / 2
  PED <- cbind(rep("FAM", nIND),
               1:nIND,
               rep(0, nIND),
               rep(0, nIND),
               rep(0, nIND),
               rep(0, nIND),
               t(Hadm))
  
  write.table(H1.table, "RASPberry_refpop1.phased", row.names = FALSE, quote = FALSE)
  write.table(H2.table, "RASPberry_refpop2.phased", row.names = FALSE, quote = FALSE)
  write.table(PED, "RASPberry.ped", col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(bim.table, "RASPberry.bim", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "      ")
  write.table(map.table, "RASPberry.map", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "      ")
}