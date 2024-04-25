######
#' @title ar.penyw
#'
#' @description Fit an autoregressive time series model to the data using penalized correlation estimation, by default selecting the order by AIC.
#'
#' @param x a univariate numeric time series.
#' @param order.max maximum lag at which to calculate the penalized/sample PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param penalized 'logical'. If 'TRUE' (the default) the penalized autoregressive model is fitted; if 'FALSE' the \code{ar} is fitted and \code{method} can be selected.
#' @param AIC 'logical', If 'TRUE' then the Akaike Information Criterion is used to choose the order of the autoregressive model. If FALSE, the model of order order.max is fitted.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param arprint  'logical'. If 'TRUE' (the default) the fitted penalized ar model is printed.
#' @param method character string specifying the method to fit the regular autoregressive time series model. Default to "yule-walker".
#' @param ... additional arguments for \code{ar} function.
#'
#' @return An object of penalized ar model fit with the following elements:
#' \describe{
#' \item{\code{order}}{The order of the fitted model. This is chosen by minimizing the AIC if 'AIC=TRUE', otherwise it is 'order.max'.}
#' \item{\code{penar}}{Estimated penalized autorregression coefficients for the fitted model.}
#' \item{\code{aic}}{The vector of AIC values.}
#' \item{\code{n.used}}{The number of observations in the time series.}
#' \item{\code{order.max}}{The value of the 'order.max' argument.}
#' \item{\code{call}}{The matched call.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' # Example for penalized ar model fit
#' ar.penyw(data)
#' }
#' @export
#####

ar.penyw=function(x, aic = TRUE, order.max = NULL, na.action = na.fail, 
                  demean = TRUE, series = NULL,...){
  if (is.null(series)) 
    series <- deparse1(substitute(x))
  ists <- is.ts(x)
  x <- na.action(as.ts(x))
  if (ists) 
    xtsp <- tsp(x)
  xfreq <- frequency(x)
  x <- as.matrix(x)
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if (any(is.na(x) != is.na(x[, 1]))) 
    stop("NAs in 'x' must be the same row-wise")
  nser <- ncol(x)
  if (demean) {
    xm <- colMeans(x, na.rm = TRUE)
    x <- sweep(x, 2L, xm, check.margin = FALSE)
  }
  else{xm <- rep.int(0, nser)}
  sampleT <- nrow(x)
  n.obs <- sum(!is.na(x[, 1]))
  order.max <- if (is.null(order.max)) 
    min(sampleT - 1L, floor(10 * log10(sampleT)))
  else floor(order.max)
  if(is.na(order.max) || order.max < 1L)
    stop("'order.max' must be >= 1")
  else if (order.max >= sampleT)
    stop("'order.max' must be less than 'n'")
  

  if(!aic){
    pencoef <- DLpencoef(x, lag.max = order.max, ...)
    penorder <- order.max
    AICpen <- NULL
  }
  else{
    penpacf=pacf(x,lag.max=order.max,demean=FALSE,plot=FALSE)
    
    AICpen <- apply(matrix(1:order.max,ncol=1), MARGIN=1,FUN=function(i){
      sampleT*log(apply(x,MARGIN=2,FUN=var)* # nser length of variances
                    apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(j){
                      prod(1-penpacf$acf[1:i,j,j]^2) # nser length of products
                    })) + 2*i
    })
    if(nser==1){
      AICpen=matrix(AICpen,nrow=1)
    }
    AICpen <- cbind(sampleT * log(apply(x,MARGIN=2,FUN=var)),AICpen)
    
    penorder <- apply(AICpen,MARGIN=1,FUN=which.min)-1
    pencoef <- lapply(1:nser,FUN=function(i){
      return(DLpencoef(x[, i], lag.max = penorder[i], ...))
    })
  }
  
  var.pred=apply(x,MARGIN=2,FUN=var)* # nser length of variances
    apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(j){
      prod(1-penpacf$acf[1:penorder[j],j,j]^2) # nser length of products
    })
  resid <- apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(j){
    x[(penorder[j]+1):sampleT,j] - 
      t(pencoef[[1]])%*%t(apply(matrix(1:penorder[j],ncol=1),MARGIN=1,FUN=function(i){x[(penorder[j]-i+1):(sampleT-i),j]}))
  })
  
  res <- list(order = penorder, ar = pencoef, var.pred = var.pred, 
              x.mean = drop(xm), aic = AICpen, n.used = sampleT, n.obs = n.obs, 
              order.max = order.max, partialacf = penpacf, resid = resid, 
              method = "Penalized Yule-Walker", series = series, frequency = xfreq, 
              call = match.call())
  class(res) <- "ar"
  return(res)
}