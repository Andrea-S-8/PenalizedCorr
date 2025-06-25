######
#' @title Auto- and Cross- Covariance and -Correlation Function Estimation
#'
#' @description It is well known that the default acf/pacf sample estimates are biased.
#' This function provides an updated version of the \code{stats::acf} function with
#' default values to calculate the penalized acf/pacf.
#' The function acf computes (and by default plots) estimates of the autocovariance 
#' or autocorrelation function. Function pacf is the function used for the partial 
#' autocorrelations. Function ccf computes the cross-correlation or cross-covariance 
#' of two univariate series.
#'
#' @param x a univariate or multivariate numeric time series object or a numeric vector or matrix.
#' @param lag.max maximum lag at which to calculate the ACF/PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param type character string giving the type of acf to be computed. Allowed values are "correlation" (the default), "covariance" or "partial". Will be partially matched.
#' @param plot  logical. If TRUE (the default) the ACF/PACF is plotted.
#' @param na.action function to be called to handle missing values. \code{na.pass} can be used.
#' @param demean logical. Should a mean be estimated and subtracted before correlations are calculated?
#' @param penalized logical. If \code{TRUE} (the default) the penalized ACF/PACF is computed; if \code{FALSE} the sample ACF/PACF is computed using \code{stats:acf}.
#' @param lh sequence of threshold values across h, default is \code{NULL}. Could be a single value (repeated for all h), a single vector of length h (repeated for all nser), or a h x nser matrix. Default is data driven choice.
#' @param estimate character vector of the estimation method for the ACF, options are \code{"direct"} (default) or \code{"invertpacf"}.  \code{"invertpacf"} is preferred when the data can be approximated by a low order AR model.
#' @param ... additional arguments for specific methods or plotting.
#'
#' @details
#' For \code{type = "correlation"} and \code{"covariance"}, if \code{penalized=FALSE} the estimates are based on the sample covariance as in \code{stats::aacf}. (The lag 0 autocorrelation is fixed at 1 by convention.)  It is well known that the 
#' sample acf/pacf estimates are biased, using \code{penalized=TRUE} results in an unbiased estimate based on shrinkage towards a target rather than uniformly shrinkage all lags towards zero.  See references for full technical details.
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
#' @return An object of class \code{"acf"} with the following elements:
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
#' ### AR(1)
#' set.seed(1234)
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' # Examples for penalized ACF/PACF and sample ACF/PACF
#' acf(data) # penalized estimate
#' acf(data, penalized = FALSE) # usual stats::acf() implementation
#' acf(x,estimate="invertpacf") # estimate the acf by inverting the pacf
#' acf(data, type ="partial") # penalized partial correlation estimate
#' acf(data, type ="partial", penalized = FALSE) # usual stats::pacf() estimate
#'
#' set.seed(1234)
#' x1 <- arima.sim(n=100, model=list(ar=0.5))
#' x2 <- arima.sim(n=100, model=list(ar=0.1))
#' x3 <- arima.sim(n=100, model=list(ar=0.9))
#' x <- cbind(x1, x2, x3)
#'
#' acf(x) #penalized estimate
#' acf(x, penalized = FALSE) # usual stats::acf() implementation
#' acf(x,estimate="invertpacf") # estimate the acf by inverting the pacf
#' acf(x, type ="partial") # penalized partial correlation estimate
#' acf(x, type ="partial", penalized = FALSE) # usual stats::pacf() estimate
#' 
#' 
#' ### MA(1)
#' set.seed(1234)
#' data <- arima.sim(n=100, model=list(ma=0.7))
#'
#' # Examples for penalized ACF/PACF and sample ACF/PACF
#' acf(data) # penalized estimate
#' acf(data, penalized = FALSE) # usual stats::acf() implementation
#' acf(x,estimate="invertpacf") # estimate the acf by inverting the pacf
#' acf(data, type ="partial") # penalized partial correlation estimate
#' acf(data, type ="partial", penalized = FALSE) # usual stats::pacf() estimate
#'
#' set.seed(1234)
#' x1 <- arima.sim(n=100, model=list(ma=0.5))
#' x2 <- arima.sim(n=100, model=list(ma=0.1))
#' x3 <- arima.sim(n=100, model=list(ma=0.9))
#' x <- cbind(x1, x2, x3)
#'
#' acf(x) #penalized estimate
#' acf(x, penalized = FALSE) # usual stats::acf() implementation
#' acf(x,estimate="invertpacf") # estimate the acf by inverting the pacf
#' acf(x, type ="partial") # penalized partial correlation estimate
#' acf(x, type ="partial", penalized = FALSE) # usual stats::pacf() estimate
#' }
#' @export
#' @importFrom stats as.ts
#' @importFrom stats na.fail
#' @importFrom stats var
#####

acf <-
function (x, lag.max = NULL, type = c("correlation", "covariance", 
    "partial"), plot = TRUE, na.action = na.fail, demean = TRUE, 
    penalized=TRUE,lh=NULL,estimate="direct",...){

    type <- match.arg(type)
    x <- na.action(as.ts(x))
    x <- as.matrix(x)
    if (!is.numeric(x)){stop("'x' must be numeric")}
    if(!is.logical(penalized)){stop("penalized must be logical")}
    sampleT <- as.integer(nrow(x))
    nser <- as.integer(ncol(x))
    
    if (type == "partial") {
        m <- match.call()
        m[[1L]] <- quote(pacf)
        m$type <- NULL
        return(eval(m, parent.frame()))
    }
    else if(type=="covariance"){
      acf=acf(x,lag.max=lag.max,type="correlation",plot=plot,na.action=na.action,
              demean=demean,penalized=penalized,lh=lh,estimate=estimate,...)
      for(i in 1:nser){
        acf$acf[,i,i]=acf$acf[,i,i]*var(x[,i])*(sampleT-1)/sampleT
      }
      return(acf)
    }
    if(estimate=="direct"){
      if(!penalized){ # not penalised so run usual acf
        acf=stats::acf(x,lag.max,type,plot=FALSE,na.action,demean,...)
        acf$penalized=FALSE
        acf$lh=NULL
      }
      else{ #run penalised estimation
        acf=corrected(x,lag.max,type,na.action,demean,lh,...)
        acf$penalized=TRUE
      }
      acf$estimate="direct"
    }
    else if(estimate=="invertpacf"){
      acf=invertpacf(x=x,lag.max=lag.max,type=type,na.action=na.action,
                     demean=demean,penalized=penalized,lh=lh,...)
      acf$penalized=penalized
      acf$estimate="invertpacf"
     }
    else{stop("The estimate argument can only take values 'direct' or 'invertpacf'.")}

    if(plot){
      extra.args=list(...)
      if(any(names(extra.args)=="ylab")){
        plot(acf,...)
      }
      else if(acf$penalized==TRUE){
        plot(acf,ylab="Penalized ACF",...)
      }
      else{
        plot(acf,...)
      }
      invisible(acf)
    }
    else{return(acf)}
}
