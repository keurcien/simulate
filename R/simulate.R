#' Simulation tools
#'
#' \code{simu.write.ancstrl} creates haplotype matrices and genotype matrices 
#' from Beagle outputs. 
#'
#' @param obj.vcfR an object of class vcfR obtained with a Beagle output.
#' @param pop a string vector specifying the population labels for 
#' each individual present in the VCF file.
#' @param ancestral.1 a character string or an integer specifying the first 
#' ancestral population.
#' @param ancestral.2 a character string or an integer specifying the second 
#' ancestral population.
#' @param recombinationRate a numerical value expressed in centiMorgans per 
#' Megabase specifying the mean recombination rate of the species.
#' @param keep a list specifying the markers to keep. If missing,
#' all markers will be kept.
#' @param output.name a character string specifying the name of the output map
#' file.
#' 
#' @importFrom vcfR read.vcfR
#' 
#' @return 
#'
#' @export
#'
simu.write.ancstrl = function(obj.vcfR, 
                              pop, 
                              ancestral.1, 
                              ancestral.2, 
                              recombinationRate, 
                              keep, 
                              output.name = "output"){
  
  hap <- extract.haps(obj.vcfR)
  
  hap.pop <- vector(length = 2 * length(pop), mode = "numeric")
  hap.pop[seq(1, length(hap.pop), by = 2)] <- pop
  hap.pop[seq(2, length(hap.pop), by = 2)] <- pop
  pos <- getPOS(obj.vcfR)
  gen <- recombinationRate * pos * 1e-6 # centiMorgans per Megabase
  
  ref <- getREF(obj.vcfR)
  alt <- getALT(obj.vcfR)
  single.alt <- which(alt %in% unique(ref))
  
  hap <- hap[single.alt, ]
  ref <- ref[single.alt]
  alt <- alt[single.alt]
  pos <- pos[single.alt]
  gen <- gen[single.alt]
  
  if (!missing(keep)){
    hap <- hap[keep, ]
    ref <- ref[keep]
    alt <- alt[keep]
    pos <- pos[keep]
    gen <- gen[keep]
  }
  
  nSNP <- nrow(hap)
  nIND <- ncol(hap) / 2
  hap.1 <- hap[, (hap.pop == ancestral.1)]
  hap.2 <- hap[, (hap.pop == ancestral.2)]
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
  
  rates <- cbind(pos, gen)
  
  cat(paste0("Writing ", output.name, ".map..."))
  write.table(rates, 
              paste0(output.name, ".map"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  cat("DONE\n")
  cat(paste0("Writing ", output.name, "_refpop1.phased..."))
  write.table(hap.1, 
              paste0(output.name, "_refpop1.phased"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  cat("DONE\n")
  cat(paste0("Writing ", output.name, "_refpop2.phased..."))
  write.table(hap.2, 
              paste0(output.name, "_refpop2.phased"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  cat("DONE\n")
  cat(paste0("Writing ", output.name, "_refpop1_integer.phased..."))
  write.table(hap.int.1, 
              paste0(output.name, "_refpop1_integer.phased"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  cat("DONE\n")
  cat(paste0("Writing ", output.name, "_refpop2_integer.phased..."))
  write.table(hap.int.2, 
              paste0(output.name, "_refpop2_integer.phased"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  cat("DONE\n")
  cat(paste0("Writing ", output.name, "_refpop1.pcadapt..."))
  write.table(geno.int.1, 
              paste0(output.name, "_refpop1.pcadapt"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  cat("DONE\n")
  cat(paste0("Writing ", output.name, "_refpop2.pcadapt..."))
  write.table(geno.int.2, 
              paste0(output.name, "_refpop2.pcadapt"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote = FALSE)
  cat("DONE\n")
}

#' Simulation tools
#'
#' \code{generate_hybrid} generates one hybrid individual using the 
#' recombination map and ancestral haplotypes.
#'
#' @param refpop1_integer a haplotype matrix. 
#' @param refpop2_integer a haplotype matrix.
#' @param refpop1 a haplotype matrix.
#' @param refpop2 a haplotype matrix.
#' @param alpha a numerical value.
#' @param beta a numerical value.
#' @param jumps a numerical vector. 
#' @param ancestry.switch a numerical vector.
#' 
#' @return 
#' 
#' @useDynLib simulate
#'
#' @export
#'
generate_one_hybrid = function(refpop1_integer, 
                               refpop2_integer,
                               refpop1,
                               refpop2,
                               alpha, 
                               beta = 1 - alpha,
                               jumps, 
                               ancestry.switch){
  mu <- alpha
  sig <- 0.1
  b_alpha <- mu * mu * ((1 - mu) / (sig * sig) - 1 / mu)
  b_beta <- b_alpha * (1 / mu - 1)
  nHAP <- min(ncol(refpop1_integer), ncol(refpop2_integer))
  if (nrow(refpop1_integer) != nrow(refpop2_integer)){
    stop("Ancestral populations should contain the same number of markers.")
  }
  nSNP <- nrow(refpop1_integer)
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
  haplotype.1.integer <- vector(mode = "numeric", length = nSNP)
  haplotype.2.integer <- vector(mode = "numeric", length = nSNP)
  haplotype.1 <- vector(mode = "character", length = nSNP)
  haplotype.2 <- vector(mode = "character", length = nSNP)
  true.ancestry <- vector(mode = "numeric", length = nSNP)
  for (i in 1:n.chunks){
    #p <- alpha
    p <- rbeta(n = 1, shape1 = b_alpha, shape2 = b_beta)
    chunk <- beg[i]:end[i]
    inter <- intersect(chunk, ancestry.switch)
    if (length(inter) > 0){
      p <- beta
    }
    nbino <- rbinom(1, 2, prob = p)
    if (nbino == 2){
      haplotype.1.integer[chunk] <- refpop1_integer[chunk, idx.father[i]]
      haplotype.2.integer[chunk] <- refpop1_integer[chunk, idx.mother[i]]
      haplotype.1[chunk] <- refpop1[chunk, idx.father[i]]
      haplotype.2[chunk] <- refpop1[chunk, idx.mother[i]]
      true.ancestry[chunk] <- 11
    } else if (nbino == 1){
      haplotype.1.integer[chunk] <- refpop1_integer[chunk, idx.father[i]]
      haplotype.2.integer[chunk] <- refpop2_integer[chunk, idx.mother[i]]
      haplotype.1[chunk] <- refpop1[chunk, idx.father[i]]
      haplotype.2[chunk] <- refpop2[chunk, idx.mother[i]]
      true.ancestry[chunk] <- 12
    } else if (nbino == 0){
      haplotype.1.integer[chunk] <- refpop2_integer[chunk, idx.father[i]]
      haplotype.2.integer[chunk] <- refpop2_integer[chunk, idx.mother[i]]
      haplotype.1[chunk] <- refpop2[chunk, idx.father[i]]
      haplotype.2[chunk] <- refpop2[chunk, idx.mother[i]]
      true.ancestry[chunk] <- 22
    }
  }
  return(list(h1.int = haplotype.1.integer, 
              h2.int = haplotype.2.integer,
              h1 = haplotype.1,
              h2 = haplotype.2,
              true.ancestry = true.ancestry))
}

#' Simulation tools
#'
#' \code{generate_hybrid_matrix} generates hybrid individuals using the 
#' recombination map and ancestral haplotypes.
#'
#' @param refpop1_integer a haplotype matrix. 
#' @param refpop2_integer a haplotype matrix.
#' @param refpop1 a haplotype matrix.
#' @param refpop2 a haplotype matrix.
#' @param alpha a numerical value.
#' @param beta a numerical value.
#' @param gen_map a numerical vector.
#' @param n.hyb an integer.
#' @param lambda a numerical value.
#' @param jumps a numerical vector. 
#' @param ancestry.switch a numerical vector.
#' 
#' @return 
#'
#' @export
#'
generate_hybrid_matrix = function(refpop1_integer,
                                  refpop2_integer,
                                  refpop1,
                                  refpop2,
                                  alpha = 0.5, 
                                  beta = 1 - alpha, 
                                  gen_map, 
                                  n.hyb, 
                                  lambda = 1.0,
                                  ancestry.switch = NULL){
  H.integer <- matrix(0, nrow = nrow(refpop1_integer), ncol = (2 * n.hyb))
  H <- matrix(0, nrow = nrow(refpop1), ncol = (2 * n.hyb))
  true.ancestry.matrix <- matrix(0, nrow = nrow(refpop1_integer), ncol = n.hyb)
  for (i in 1:n.hyb){
    jumps <- jumps_from_map(gen_map * 1e-2, lambda = lambda) # g in Morgans
    h <- generate_one_hybrid(refpop1_integer, 
                             refpop2_integer, 
                             refpop1,
                             refpop2,
                             alpha, 
                             beta, 
                             jumps, 
                             ancestry.switch = ancestry.switch)
    H.integer[, (2 * i - 1)] <- h$h1.int
    H.integer[, (2 * i)] <- h$h2.int
    H[, (2 * i - 1)] <- h$h1
    H[, (2 * i)] <- h$h2
    true.ancestry.matrix[, i] <- h$true.ancestry
  }
  return(list(H.integer = H.integer, 
              H = H, 
              true.ancestry.matrix = true.ancestry.matrix))
}

#' Simulation tools
#'
#' \code{check_intersect}
#'
#' @param p a list containing two integers.
#' @param ground.truth a vector of integers.
#' 
#' @return 
#'
#' @export
#'
check_intersect = function(p, ground.truth){
  for (i in length(p$beg)){
    if (length(intersect(p$beg[i]:p$end[i], ground.truth)) > 0){
      return(TRUE)    
    }  
  }
  return(FALSE)
}

#' Simulation tools
#'
#' \code{get_dist_pop}
#'
#' @param H1 a haplotype matrix. 
#' @param H2 a haplotype matrix.
#' 
#' @return 
#'
#' @export
#'
get_dist_pop = function(H1, H2){
  nSNP <- nrow(H1)
  sum.1 <- apply(H1, MARGIN = 1, sum)
  sum.2 <- apply(H2, MARGIN = 1, sum)
  dist <- abs(sum.2 - sum.1)
  return(dist)
}
