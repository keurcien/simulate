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

create_input_hapmix = function(H1, H2, H, gen_map, phys_map, chr, 
                               ancstrl.1 = "ANC1", ancstrl.2 = "ANC2", 
                               admxd = "AA", rates.digits = 8){
  cat("Creating the genotype file for ancestral population 1...")
  write.table(t(H1), paste0(ancstrl.1, ".lfmm"), col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(ancstrl.1, ".lfmm"), paste0(ancstrl.1, ".geno"), force = TRUE)
  cat("Creating the genotype file for ancestral population 2...")
  write.table(t(H2), paste0(ancstrl.2, ".lfmm"), col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(ancstrl.2, ".lfmm"), paste0(ancstrl.2, ".geno"), force = TRUE)
  cat("Creating the genotype file for the admixed population...")
  G <- haplo_to_geno(H)
  write.table(t(G), paste0(admxd, ".lfmm"), col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  LEA::lfmm2geno(paste0(admxd, ".lfmm"), paste0(admxd, ".geno"), force = TRUE)
  
  dt <- array(dim = c(nrow(H1), ncol = 4))
  dt[, 1] <- paste0("rs", 1:nrow(H1))
  dt[, 2] <- chr
  dt[, 3] <- format(round(gen_map, rates.digits), nsmall = rates.digits)
  dt[, 4] <- phys_map
  cat("Creating the snpfiles...")
  write.table(dt, paste0(ancstrl.1, "snpfile"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(dt, paste0(ancstrl.2, "snpfile"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(dt, paste0(admxd, "snpfile"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  dt.ind <- as.data.frame(matrix(nrow = ncol(G), ncol = 3))
  dt.ind[[1]] <- paste0(rep("IND_", ncol(G)), 1:ncol(G))
  dt.ind[[2]] <- "F"
  dt.ind[[3]] <- "Hybrid"
  write.table(dt.ind, paste0(admxd, ".ind"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  cat("DONE\n")
  cat("Creating the rates file...")
  l1 <- paste0(":sites:", nrow(G))
  l2 <- stringr::str_c(phys_map, collapse = " ")
  l3 <- stringr::str_c(format(gen_map, scientific = FALSE, nsmall = rates.digits), 
                       collapse = " ")
  con <- file("rates")
  writeLines(paste0(l1, "\n", l2, "\n", l3), con = con)
  close(con = con)
  cat("DONE\n")
}

create_input_eila = function(H1, H2, H, position){
  obj.eila <- NULL
  obj.eila$anc1 <- haplo_to_geno(H1)  
  obj.eila$anc2 <- haplo_to_geno(H2)
  obj.eila$admixed <- haplo_to_geno(H)
  obj.eila$position <- position
  return(obj.eila)
}

create_input_rfmix = function(H1, H2, H3, gen_map){
  G <- cbind(H1, H2, H3)
  write.table(t(G), "rfmix_alleles.lfmm", col.names = FALSE, row.names = FALSE)
  LEA::lfmm2geno("rfmix_alleles.lfmm", "rfmix_alleles.geno", force = TRUE)
  formated.gen_map <- format(gen_map * 10000, scientific = FALSE, nsmall = 8) # in cM
  write.table(formated.gen_map, "rfmix_markerlocation.txt", sep = "\n", col.names = FALSE, row.names = FALSE, quote = FALSE)
  pop <- c(rep(1, ncol(H1)), rep(2, ncol(H2)), rep(0, ncol(H3)))
  l3 <- stringr::str_c(pop, collapse = " ")
  con <- file("rfmix_classes.txt")
  writeLines(l3, con = con)
  close(con = con)
}