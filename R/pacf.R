######
#' @title pacf
#'
#' @description The Penalized Partial Autocorrelation Function (PACF) Estimation
#'
#' @param x a univariate or multivariate numeric time series.
#' @param lag.max maximum lag at which to calculate the penalized/sample PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param C a positive real number for l_h
#' @param penalized 'logical'. If 'TRUE' (the default) the penalized PACF is computed; if 'FALSE' the sample PACF is computed.
#' @param print  'logical'. If 'TRUE' (the default) the penalized/sample PACF is printed.
#' @param plot  'logical'. If 'TRUE' (the default) the penalized/sample PACF is plotted.
#' @param series names for the series. Defaults to 'deparse(substitute(x))'.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param demean 'logical'. Should a mean be estimated during estimating.
#' @param ... additional arguments for specific methods.
#'
#' @return An object of penalized/sample PACF estimation, which is a list with the following elements:
#' \describe{
#' \item{\code{acf}}{An array containing the estimated penalized/sample PACF.}
#' \item{\code{lag}}{An array containing the lags at which the PACF is estimated.}
#' \item{\code{n.used}}{The number of observation in the time series.}
#' \item{\code{series}}{The name of the time series.}
#' \item{\code{snames}}{The series names for a multivariate time series.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' # Examples for penalized PACF and sample PACF
#' pacf(data)
#' pacf(data, penalized = FALSE)
#' }
#' @export
#####

pacf <-
function (x, lag.max=NULL, plot=TRUE, na.action=na.fail, demean=TRUE,penalized=TRUE,lh=NULL,...){
  if(!penalized){ # not penalised so use standard pacf
    pacf=stats:pacf(x,lag.max,plot,na.action,...)
    pacf$penalized=FALSE
  }
  else{ # penalised output
    pacf=corrected(x,lag.max,type="partial",na.action,demean,lh,...)
    pacf$penalized=TRUE
  }
  return(pacf)
}
