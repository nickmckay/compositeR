#' simulateAutoCorrelatedUncertainty
#' @description generate synthetic timeseries
#'
#' @param sd standard deviation of the output vector
#' @param n length of output vector
#' @param mean mean of the output vector
#' @param ar Autocorrelation coefficient to use for modelling uncertainty, what fraction of the uncertainties are autocorrelated? (default = sqrt(0.5); or 50 percent autocorrelated uncertainty)
#' @param arima.order Order to use for ARIMA model used in modelling uncertainty (default = c(1,0,0))
#'
#' @return ts simulated from a from an ARIMA model with a defined mean and variance
#' @export
simulateAutoCorrelatedUncertainty <- function(sd, n, mean = 0, ar = sqrt(.5),arima.order = c(1,0,0)){
  unc <- stats::arima.sim(list(order = arima.order, ar = ar), n = n)
  unc <- (scale(unc,center = TRUE, scale = TRUE) * sd) + mean
  return(unc)
}



