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
#' \item{\code{acf}}{A max.lag x nseries x nseries array containing the estimated penalized acf/pacf.}
#' \item{\code{type}}{Character vector returning the type argument requested.}
#' \item{\code{n.used}}{Numeric of the number of points used for estimation after na.action has been applied.}
#' \item{\code{lag}}{A max.lag x nseries x nseries array containing the lags at which the acf/pacf is estimated.}
#' \item{\code{series}}{The name of the time series.}
#' \item{\code{snames}}{The series names for a multivariate time series.}
#' \item{\code{penalized}}{Logical indicating if the acf/pacf returned is penalized.}
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
#' @importFrom stats na.fail
#####

pacf <-
function (x, lag.max=NULL, plot=TRUE, na.action=na.fail, demean=TRUE,penalized=TRUE,lh=NULL,...){
  if(!penalized){ # not penalised so use standard pacf
    pacf=stats::pacf(x,lag.max,plot,na.action,...)
    pacf$penalized=FALSE
  }
  else{ # penalised output
    pacf=corrected(x,lag.max,type="partial",na.action,demean,lh,...)
    pacf$penalized=TRUE
  }
  return(pacf)
}