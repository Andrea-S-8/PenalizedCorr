######
#' @title DLpencoef
#'
#' @description Compute coefficients via the Durbin Levinson algorithm using the penalized partial autocorrelation function (PACF) estimation to obtain an penalized autocorrelation function (ACF) estimation which is positive definite.
#'
#' @param x a univariate or multivariate numeric time series.
#' @param lag.max maximum lag at which to calculate the coefficients. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param penalized 'logical'. If 'TRUE' (the default) the penalized PACF is computed; if 'FALSE' the sample PACF is computed.
#' @param lh sequence of threshold values across h. Could be a single value (repeated for all h), a single vector of length lag.max (repeated for all nser), or a h x nser matrix. Default is data driven.
#' @param return.mat 'logical'. If 'TRUE' the return form is a matrix.  If 'FALSE' the return form is an acf object.
#' @param ... additional arguments for penalized PACF estimation.
#'
#' @return A coefficients vector to estimate a penalized autocorrelation function which is positive definite.
#'
#'
#' @references Gallagher, C., Killick, R., Tan, X. (2024+) Penalized M-estimation 
#' for Autocorrelation. \emph{Submitted.}
#' 
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' DLpencoef(data)
#' DLpencoef(data, lag.max=10)
#' }
#' @keywords internal
#' @importFrom stats as.ts
#' @importFrom stats na.fail
#####

DLpencoef <- function(x, lag.max = NULL, na.action=na.fail,penalized=TRUE,lh=NULL,return.mat=FALSE,...){
  x <- na.action(as.ts(x))
  x <- as.matrix(x)
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  sampleT <- nrow(x)
  nser=ncol(x)
  
  if (is.null(lag.max)) 
    lag.max <- floor(10 * (log10(sampleT) - log10(nser)))
  
  penpacf <- pacf(x, lag.max = lag.max, plot = FALSE, na.action=na.action,penalized=penalized,lh=lh,...)
  coef <- matrix(0,nrow=lag.max, ncol=lag.max)
  coef[1,1] <- penpacf$acf[1]
  if(lag.max == 1){return(coef[1])}
  for(i in 2:lag.max){
    if(i>sqrt(sampleT)){
      penpacf <- pacf(x, lag.max = i, plot = FALSE, na.action=na.action, penalized=penalized,lh=lh[1:i,],...)
    }
    coef[i,i] <- penpacf$acf[i]
    coef[1:(i-1),i] <- coef[1:(i-1),i-1] - penpacf$acf[i] * coef[(i-1):1,i-1]
  }
  if(return.mat){
    return(coef)
  }
  else{return(coef[,ncol(coef)])}
}
