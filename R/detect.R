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
