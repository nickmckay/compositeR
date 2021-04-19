#
#' Composite with uncertainty
#'
#' @param fTS TS object of the timeseries that you want to composite
#' @param binvec A vector of bin edges for binning
#' @param spread An option to spread values over adjacent empty years (needs a better definition here)
#' @param stanFun Function to use for standardization
#' @param ageVar Name of the time vector in fTS
#' @param gaussianizeInput Optionally gaussianize input before compositing
#' @param alignInterpDirection Optionally invert timeseries with negative interpretation_direction
#' @param binFun Function to use for binning
#'
#' @return
#' @export
#'
#' @examples
compositeEnsembles <- function(fTS,
                               binvec,
                               spread = TRUE,
                               stanFun = standardizeMeanIteratively,
                               ageVar = "age",
                               gaussianizeInput = FALSE,
                               alignInterpDirection = TRUE,
                               binFun = sampleEnsembleThenBinTs,
                               ...){
  binAges <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

  binMatR <- as.matrix(purrr::map_dfc(fTS,binFun,binvec,ageVar = ageVar,spread = spread,gaussianizeInput = gaussianizeInput,alignInterpDirection = alignInterpDirection))
  binMatR[is.nan(binMatR)] <- NA

  compMat <- stanFun(ages = binAges,pdm = binMatR,...)

  #which records contributed?
  gm <- is.finite(compMat)
  comp <- rowMeans(compMat,na.rm = TRUE)
  count <- colSums(gm)
  contributed <- which(count > 0)
  return(list(composite = comp, count = count, contributed = contributed))
}






