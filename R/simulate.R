#' Simulation tools
#'
#' \code{simu.write.ancstrl} creates haplotype matrices and genotype matrices 
#' from Beagle outputs. 
#'
#' @param vcf.file a character string specifying the name of a VCF file 
#' generated with Beagle containing phased haplotypes.
#' @param pop.file a character string specifying the name of a text file 
#' containing the population labels for each individual present in the VCF file.
#' @param ancstrl.1 a character string or an integer specifying the first 
#' ancestral population.
#' @param ancstrl.2 a character string or an integer specifying the second 
#' ancestral population.
#' @param recombinationRate a numerical value expressed in cM, specifying the
#' mean recombination rate of the species.
#' @param keep an integer specifying the number of markers to keep. If missing,
#' all markers will be kept.
#' @param output.name a character string specifying the name of the output map
#' file.
#' 
#' @return 
#'
#' @export
#'
simu.write.ancstrl = function(vcf.file, pop.file, ancstrl.1, ancstrl.2, 
                              recombinationRate, keep, output.name = "output"){
  
  beagle <- read.vcfR(vcf.file)
  pop <- read.table(pop.file)[, 1]
  
  hap.pop <- vector(length = 2 * length(pop), mode = "numeric")
  hap.pop[seq(1, length(hap.pop), by = 2)] <- pop
  hap.pop[seq(2, length(hap.pop), by = 2)] <- pop
  
  pos <- getPOS(beagle)
  gen <- recombinationRate * pos / 1000000
  rates <- cbind(pos, gen)
  
  hap <- extract.haps(beagle)
  ref <- getREF(beagle)
  alt <- getALT(beagle)
  
  nSNP <- nrow(hap)
  nIND <- ncol(hap) / 2
  hap.1 <- hap[, (hap.pop == ancstrl.1)]
  hap.2 <- hap[, (hap.pop == ancstrl.2)]
  hap.int.1 <- array(0, dim = dim(hap.1))
  hap.int.2 <- array(0, dim = dim(hap.2))
  geno.int.1 <- array(0, dim = c(nrow(hap.int.1), ncol(hap.int.1) / 2))
  geno.int.2 <- array(0, dim = c(nrow(hap.int.2), ncol(hap.int.2) / 2))
  
  for (j in 1:ncol(hap.1)){
    hap.int.1[, j] <- 1 * (alt == hap.1[, j])
  }
  
  for (j in 1:ncol(hap.2)){
    hap.int.2[, j] <- 1 * (alt == hap.2[, j])
  }
  
  for (j in 1:ncol(geno.int.1)){
    geno.int.1[, j] <- hap.int.1[, (2 * j - 1)] + hap.int.1[, (2 * j)]
  }
  
  for (j in 1:ncol(geno.int.2)){
    geno.int.2[, j] <- hap.int.2[, (2 * j - 1)] + hap.int.2[, (2 * j)]
  }
  
  if (!missing(keep)){
    hap.int.1 <- hap.int.1[1:keep, ]
    hap.int.2 <- hap.int.2[1:keep, ]
    geno.int.1 <- geno.int.1[1:keep, ]
    geno.int.2 <- geno.int.2[1:keep, ]
    rates <- rates[1:keep, ]
  }
  
  write.table(rates, paste0(output.name, ".map"), col.names = FALSE, row.names = FALSE)
  cat("Writing H_1...")
  write.table(hap.int.1, "H_1", col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  cat("Writing H_2...")
  write.table(hap.int.2, "H_2", col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  cat("Writing G_1...")
  write.table(geno.int.1, "G_1", col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
  cat("Writing G_2...")
  write.table(geno.int.2, "G_2", col.names = FALSE, row.names = FALSE)
  cat("DONE\n")
}

#' Simulation tools
#'
#' \code{generate_hybrid} generates one hybrid individual using the 
#' recombination map and ancestral haplotypes.
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' @param map a numerical vector containing the genetic positions.
#' @param global.ancestry a numerical value between \code{0} and \code{1}
#' setting the average proportion of markers coming from the \code{H1}.
#' 
#' @return 
#'
#' @export
#'
generate_hybrid = function(H1, H2, map, n.gen = 10, global.ancestry = 0.5){
  nHAP <- ncol(H1)
  jumps <- jumps_builder(map = map, lambda = n.gen)
  boolean <- jumps_logical(jumps)
  p <- rbinom(1, size = 2, prob = global.ancestry)
  if (p == 2){
    idx.father <- sample(1:nHAP, size = 1)
    idx.mother <- sample(1:nHAP, size = 1)
    hap.father <- H1[, idx.father]
    hap.mother <- H1[, idx.mother]
  } else if (p == 1){
    idx.father <- sample(1:nHAP, size = 1)
    idx.mother <- sample(1:nHAP, size = 1)
    hap.father <- H1[, idx.father]
    hap.mother <- H2[, idx.mother]
  } else if (p == 0){
    idx.father <- sample(1:nHAP, size = 1)
    idx.mother <- sample(1:nHAP, size = 1)
    hap.father <- H2[, idx.father]
    hap.mother <- H2[, idx.mother]
  }
  adm.mat <- array(0, dim = c(2 * n.hyb, length(jumps)))
}

gen_one_hybrid = function(H1, H2, alpha, jumps){
  nHAP <- min(ncol(H1), ncol(H2))
  if (nrow(H1) != nrow(H2)){
    stop("Ancestral populations should contain the same number of markers.")
  }
  nSNP <- nrow(H1)
  n.jumps <- sum(jumps)
  n.chunks <- n.jumps + 1 
  beg <- vector(mode = "numeric", length = n.chunks)
  end <- vector(mode = "numeric", length = n.chunks)
  beg[1] <- 1
  end[n.chunks] <- nSNP
  jumps.loc <- which(jumps == 1)
  beg[-1] <- jumps.loc
  end[1:n.jumps] <- pmax(jumps.loc - 1, 1)
  idx.father <- sample(1:nHAP, size = n.chunks, replace = TRUE)
  idx.mother <- sample(1:nHAP, size = n.chunks, replace = TRUE)
  haplotype.1 <- vector(mode = "numeric", length = nSNP)
  haplotype.2 <- vector(mode = "numeric", length = nSNP)
  for (i in 1:n.chunks){
    p <- 1 - alpha
    nbino <- rbinom(1, 2, prob = p)
    if (nbino == 2){
      haplotype.1[beg[i]:end[i]] <- H1[beg[i]:end[i], idx.father[i]]
      haplotype.2[beg[i]:end[i]] <- H1[beg[i]:end[i], idx.mother[i]]
    } else if (nbino == 1){
      haplotype.1[beg[i]:end[i]] <- H1[beg[i]:end[i], idx.father[i]]
      haplotype.2[beg[i]:end[i]] <- H2[beg[i]:end[i], idx.mother[i]]
    } else if (nbino == 0){
      haplotype.1[beg[i]:end[i]] <- H2[beg[i]:end[i], idx.father[i]]
      haplotype.2[beg[i]:end[i]] <- H2[beg[i]:end[i], idx.mother[i]]
    }
  }
  return(list(h1 = haplotype.1, h2 = haplotype.2))
}

gen_hybrid_matrix = function(H1, H2, alpha = 0.5, gen_map, n.hyb = ncol(H1) / 2, lambda = 1.0){
  H <- matrix(0, nrow = nrow(H1), ncol = (2 * n.hyb))
  for (i in 1:n.hyb){
    jumps <- jumps_from_map(gen_map, lambda = lambda)  
    h <- gen_one_hybrid(H1, H2, alpha, jumps)
    H[, (2 * i - 1)] <- h$h1
    H[, (2 * i)] <- h$h2
  }
  return(H)
}
