#composite with uncertainty
compositeEnsembles <- function(fTS,binvec,spread = TRUE,stanFun = standardizeMeanIteratively,ageVar = "age",gaussianizeInput = FALSE,alignInterpDirection = TRUE,binFun = sampleEnsembleThenBinTs,...){
  binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

  binMatR <- as.matrix(purrr::map_dfc(fTS,binFun,binvec,ageVar = ageVar,spread = spread,gaussianizeInput = gaussianizeInput,alignInterpDirection = alignInterpDirection))
  binMatR[is.nan(binMatR)] <- NA

  compMat <- stanFun(ages = binAges,pdm = binMatR,...)

  #which records contributed?
  gm <- is.finite(compMat)
  comp <- rowMeans(compMat,na.rm = TRUE)
  count <- rowSums(gm)
  contributed <- which(count > 0)
  return(list(composite = comp, count = count, contributed = contributed))
}






