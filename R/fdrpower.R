#' @export
#'
compute.fdr = function(list, ground.truth){
  if (length(list) == 0){
    x <- 0
  } else {
    x <- sum(!(list %in% ground.truth)) / length(list)
  }
  return(x)
}

#' @export
#'
compute.power = function(list, ground.truth){
  if (length(ground.truth) == 0){
    warning("The list of true positives is empty.")
  } else {
    x <- sum(list %in% ground.truth) / length(ground.truth)
  }
  return(x)
}

#' @export
#'
create.fdr.pow = function(list, ground.truth, lmax = length(list), smooth = TRUE){
  fdr <- 0
  pow <- 0
  s <- seq(1, lmax, by = 1)
  for (k in s){
    l <- list[1:k]
    fdr <- c(fdr, compute.fdr(l, ground.truth))
    pow <- c(pow, compute.power(l, ground.truth))
  }
  res <- NULL
  idx <- sort(fdr, decreasing = FALSE, index.return = TRUE)$ix
  fdr <- fdr[idx]
  pow <- pow[idx]
  if (smooth == FALSE){
    res$fdr <- fdr
    res$pow <- pow
  } else if (smooth == TRUE){
    pow.old <- 0
    fdr.old <- 0
    aux.fdr <- 0
    aux.pow <- 0
    for (val in 1:length(fdr)){
      if (pow[val] >= pow.old){
        aux.fdr <- c(aux.fdr, fdr[val])
        aux.pow <- c(aux.pow, pow[val])
        pow.old <- pow[val]
      }
    }
    res$fdr <- aux.fdr
    res$pow <- aux.pow
  }
  return(res)
}

#' @export
#'
create.df = function(res, soft.name){
  ncol <- length(res$fdr)
  df <- array(0, dim = c(ncol, 3))
  df[, 1] <- rep(soft.name, ncol)
  df[, 2] <- res$fdr
  df[, 3] <- res$pow
  return(as.data.frame(df))
}