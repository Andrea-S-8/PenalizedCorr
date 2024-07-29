######
#' @title pacf
#'
#' @description The Penalized Partial Autocorrelation Function (PACF) Estimation
#'
#' @param x a univariate or multivariate numeric time series.
#' @param lag.max maximum lag at which to calculate the penalized/sample PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param plot  'logical'. If 'TRUE' (the default) the penalized/sample PACF is plotted.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param demean 'logical'. Should a mean be estimated during estimating.
#' @param penalized 'logical'. If 'TRUE' (the default) the penalized PACF is computed; if 'FALSE' the sample PACF is computed.
#' @param lh sequence of threshold values across h. Could be a single value (repeated for all h), a single vector of length h (repeated for all nser), or a h x nser matrix. Default is data driven.
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
#' # AR(1)
#' set.seed(1234)
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' pacf(data) # penalized estimate
#' pacf(data, penalized = FALSE) # standard stats::pacf() estimate
#' 
#' set.seed(1234)
#' x1 <- arima.sim(n=100, model=list(ma=0.5))
#' x2 <- arima.sim(n=100, model=list(ma=0.1))
#' x3 <- arima.sim(n=100, model=list(ma=0.9))
#' x <- cbind(x1, x2, x3)
#' 
#' pacf(x) # penalized estimate
#' pacf(x, penalized = FALSE) # standard stats::pacf() estimate
#' 
#' # MA(1)
#' ### MA(1)
#' set.seed(1234)
#' data <- arima.sim(n=100, model=list(ma=0.7))
#' 
#' pacf(data) # penalized pacf
#' pacf(data, penalized = FALSE) # stats::pacf
#' 
#' set.seed(1234)
#' x1 <- arima.sim(n=100, model=list(ma=0.5))
#' x2 <- arima.sim(n=100, model=list(ma=0.1))
#' x3 <- arima.sim(n=100, model=list(ma=0.9))
#' x <- cbind(x1, x2, x3)
#' 
#' pacf(x) # penalized pacf
#' pacf(x, penalized = FALSE) # stats::pacf
#' 
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
  if(plot){
    plot(pacf, ...)
    invisible(pacf)
  }
  else return(pacf)
}