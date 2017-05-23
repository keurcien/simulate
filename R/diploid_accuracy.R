getGenotype = function(row){
  row = row[order(row)]
  return(paste(row, collapse = ""))
}

rfmix.diploid.accuracy <- function(gt, rfmix.output){
  #solution = read.table(gt, stringsAsFactors = F, header = F)
  solution = gt
  numAdm = ncol(solution)
  output = read.table(rfmix.output, stringsAsFactors = F)
  numSnps = nrow(solution)
  totalLocs = numAdm * numSnps
  newTotal = 0
  for(adm in 1:numAdm){
    called = output[, c(1 + 2 * (adm - 1), 2 * (adm - 1) + 2)]
    true = solution[, adm]
    calledGen = apply(called, 1, getGenotype)
    newAccuracy = length(which(calledGen == as.character(true)))
    newTotal = newTotal + newAccuracy
    cat(newAccuracy / numSnps, "\n")
  }
  cat("Mean Accuracy", newTotal / totalLocs, "\n")
}

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






