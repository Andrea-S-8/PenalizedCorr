######
#' @title ar.penyw
#'
#' @description Fit an autoregressive time series model to the data using penalized correlation estimation, by default selecting the order by AIC.
#'
#' @param x a univariate numeric time series.
#' @param aic 'logical', If 'TRUE' then the Akaike Information Criterion is used to choose the order of the autoregressive model. If FALSE, the model of order order.max is fitted.
#' @param order.max maximum lag at which to calculate the penalized/sample PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param demean 'logical'. Should a mean be estimated during estimating.
#' @param series names for the series. Defaults to \code{deparse1(substitute(x))}.
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
#' @importFrom stats na.fail
#' @importFrom stats is.ts
#' @importFrom stats as.ts
#' @importFrom stats tsp
#' @importFrom stats frequency
#' @importFrom stats var
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
      if(penorder[i]==0){
        return(NULL)
      }
      return(DLpencoef(x[, i], lag.max = penorder[i], ...))
    })
  }
  
  var.pred=apply(x,MARGIN=2,FUN=var)* # nser length of variances
    apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(j){
      if(penorder[j]==0){return(1)}
      prod(1-penpacf$acf[1:penorder[j],j,j]^2) # nser length of products
    })
  resid <- apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(j){
    if(penorder[j]==0){return(x[,j])} # independent so resid=data
    x[(penorder[j]+1):sampleT,j] - 
      t(pencoef[[j]])%*%t(apply(matrix(1:penorder[j],ncol=1),MARGIN=1,FUN=function(i){x[(penorder[j]-i+1):(sampleT-i),j]}))
  })
  if(nser==1){pencoef=pencoef[[1]]} # so the print for AR works.
  
  res <- list(order = penorder, ar = pencoef, var.pred = var.pred, 
              x.mean = drop(xm), aic = AICpen, n.used = sampleT, n.obs = n.obs, 
              order.max = order.max, partialacf = penpacf, resid = resid, 
              method = "Penalized Yule-Walker", series = series, frequency = xfreq, 
              call = match.call())
  class(res) <- "ar"
  return(res)
}