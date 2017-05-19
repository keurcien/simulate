detect_outlier_region = function(stat, n.sd){
  s <- sd(stat[!is.na(stat)])
  idx <- which(stat > (n.sd * s))
  list.beg.reg <- idx[1]
  list.end.reg <- NULL
  for (i in 2:length(idx)){
    if ((idx[i] - idx[i - 1]) > 1){
      list.end.reg <- c(list.end.reg, idx[i - 1])
      list.beg.reg <- c(list.beg.reg, idx[i])
    }
  }
  list.end.reg <- c(list.end.reg, tail(idx, n = 1))
  return(list(beg = list.beg.reg, end = list.end.reg))
}

check_intersect = function(p, ground.truth){
  for (i in length(p$beg)){
    if (length(intersect(p$beg[i]:p$end[i], ground.truth)) > 0){
      return(TRUE)    
    }  
  }
  return(FALSE)
}

haplo_to_ancestry = function(H, anc){
  nSNP <- nrow(H)
  nIND <- ncol(H)
  mean_ancestry <- vector(mode = "numeric", length = nrow(H))
  
  if (anc == 1){
    mean_ancestry <- mean_ancestry + apply(solution, MARGIN = 1, FUN = function(x){sum(x == "12" || x == "21")})   
    mean_ancestry <- mean_ancestry + 2 * apply(solution, MARGIN = 1, FUN = function(x){sum(x == "11")})   
  } else if (anc == 2){
    mean_ancestry <- mean_ancestry + apply(solution, MARGIN = 1, FUN = function(x){sum(x == "12" || x == "21")})
    mean_ancestry <- mean_ancestry + 2 * apply(solution, MARGIN = 1, FUN = function(x){sum(x == "22")})   
  }
  
  return(mean_ancestry / (2 * nIND))
}