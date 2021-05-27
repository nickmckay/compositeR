#' simulateAutoCorrelatedUncertainty
#'
#' @param sd
#' @param n
#' @param mean
#' @param ar
#' @param arima.order
#'
#' @return
#' @export
#'
#' @examples
simulateAutoCorrelatedUncertainty <- function(sd, n, mean = 0, ar = sqrt(.5),arima.order = c(1,0,0)){
  unc <- arima.sim(list(order = arima.order, ar = ar), n = n)
  unc <- (scale(unc,center = TRUE, scale = TRUE) * sd) + mean
  return(unc)
}



