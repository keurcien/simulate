#' Ancestry tools
#'
#' \code{haplo_to_ancestry} computes a vector of local ancestries from a 
#' haplotype matrix.
#'
#' @param H1 a haplotype matrix. 
#' @param anc an integer.
#' 
#' @return 
#'
#' @export
#'
haplo_to_ancestry = function(H, anc){
  nSNP <- nrow(H)
  nIND <- ncol(H)
  mean_ancestry <- vector(mode = "numeric", length = nrow(H))
  if (anc == 1){
    mean_ancestry <- mean_ancestry + apply(H, MARGIN = 1, FUN = function(x){sum(x %in% c("12", "21"))})   
    mean_ancestry <- mean_ancestry + 2 * apply(H, MARGIN = 1, FUN = function(x){sum(x == "11")})   
  } else if (anc == 2){
    mean_ancestry <- mean_ancestry + apply(H, MARGIN = 1, FUN = function(x){sum(x %in% c("12", "21"))})
    mean_ancestry <- mean_ancestry + 2 * apply(H, MARGIN = 1, FUN = function(x){sum(x == "22")})   
  }
  return(mean_ancestry / (2 * nIND))
}

#' Ancestry tools
#'
#' \code{rfmix.local.ancestry} computes a vector of local ancestries from an
#' output file obtained with RFMix.
#'
#' @param rfmix.output a character string. 
#' 
#' @return 
#'
#' @export
#'
rfmix.local.ancestry <- function(rfmix.output){
  output <- read.table(rfmix.output, stringsAsFactors = F)
  numAdm <- ncol(output) / 2
  numSnps <- nrow(output)
  mat <- matrix(0, numSnps, numAdm)
  for(adm in 1:numAdm){
    mat[, adm] = paste0(output[, 2 * adm - 1], output[, 2 * adm])
  }
  return(mat)
}

#' Ancestry tools
#'
#' \code{display.ancestry}
#'
#' @param rfmix.output a character string. 
#' 
#' @return 
#'
#' @export
#'
display.ancestry <- function(ancestry.matrix){
    im <- matrix(0, nrow(ancestry.matrix), ncol(ancestry.matrix))
    im[ancestry.matrix == "22"] <- 2
    im[ancestry.matrix == "11"] <- 0
    im[ancestry.matrix == "12"] <- 1
    im[ancestry.matrix == "21"] <- 1
    return(im)
    #image(im, axes = FALSE, xlab = "", ylab = "")
}



