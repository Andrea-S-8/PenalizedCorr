######
#' @title Partial Auto- and Cross- Covariance and -Correlation Function Estimation
#'
#' @description It is well known that the default pacf sample estimates are biased.
#' This function provides an updated version of the \code{stats::pacf} function with
#' default values to calculate the penalized pacf.
#' The function pacf computes (and by default plots) estimates of the partial 
#' autocovariance or autocorrelation function.
#' 
#' @param x a univariate or multivariate numeric time series object or a numeric vector or matrix.
#' @param lag.max maximum lag at which to calculate the ACF/PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param plot  logical. If TRUE (the default) the ACF/PACF is plotted.
#' @param na.action function to be called to handle missing values. \code{na.pass} can be used.
#' @param demean logical. Should a mean be estimated and subtracted before correlations are calculated?
#' @param penalized logical. If \code{TRUE} (the default) the penalized ACF/PACF is computed; if \code{FALSE} the sample ACF/PACF is computed using \code{stats:acf}.
#' @param lh sequence of threshold values across h, default is \code{NULL}. Could be a single value (repeated for all h), a single vector of length h (repeated for all nser), or a h x nser matrix. Default is data driven choice.
#' @param ... additional arguments for specific methods or plotting.
#'
#' @details
#' For \code{type = "correlation"} and \code{"covariance"}, if \code{penalized=FALSE} the estimates are based on the sample covariance as in \code{stats::pacf}. It is well known that the 
#' sample pacf estimates are biased, using \code{penalized=TRUE} results in an unbiased estimate based on shrinkage towards a target rather than uniformly shrinkage all lags towards zero.  See references for full technical details.
#' 
#' By default, no missing values are allowed. If the \code{na.action} function passes through missing values (as \code{na.pass} does), the covariances are computed from the complete cases. This means that the estimate computed may well not 
#' be a valid autocorrelation sequence, and may contain missing values. Missing values are not allowed when computing the PACF of a multivariate time series.
#' 
#' The partial correlation coefficient is estimated by fitting autoregressive models of successively higher orders up to \code{lag.max}.
#' 
#' The generic function \code{plot} has a method for objects of class \code{"acf"}.
#' 
#' The lag is returned and plotted in units of time, and not numbers of observations.
#' 
#' There are \code{print} and subsetting methods for objects of class \code{"acf"}.
#' 
#'@return An object of class \code{"acf"} with the following elements:
#' \describe{
#' \item{\code{acf}}{A \code{lag.max} x \code{nseries} x \code{nseries} array containing the estimated penalized acf/pacf.}
#' \item{\code{type}}{Character vector returning the \code{type} argument.}
#' \item{\code{n.used}}{Numeric of the number of points used for estimation after \code{na.action} has been applied.}
#' \item{\code{lag}}{A \code{lag.max} x \code{nseries} x \code{nseries} array containing the lags at which the acf/pacf is estimated.}
#' \item{\code{series}}{The name of the time series, \code{x}.}
#' \item{\code{snames}}{The series names for a multivariate time series.}
#' \item{\code{penalized}}{Logical returning the \code{penalized} argument.}
#' \item{\code{estimate}}{Character vector returning the \code{estimate} argument.}
#' }
#' 
#' @references Gallagher, C., Killick, R., Tan, X. (2024+) Penalized M-estimation 
#' for Autocorrelation. \emph{Submitted.}
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
  if(!is.logical(penalized)){stop("penalized must be logical")}
  if(!penalized){ # not penalised so use standard pacf
    pacf=stats::pacf(x,lag.max,plot,na.action,...)
    pacf$penalized=FALSE
  }
  else{ # penalised output
    pacf=corrected(x,lag.max,type="partial",na.action,demean,lh,...)
    pacf$penalized=TRUE
  }
  if(plot){
    extra.args=list(...)
    if(any(names(extra.args)=="ylab")){
      plot(pacf,...)
    }
    else if(pacf$penalized==TRUE){
      plot(pacf,ylab="Penalized PACF",...)
    }
    else{
      plot(pacf,...)
    }
    invisible(pacf)
  }
  else return(pacf)
}
