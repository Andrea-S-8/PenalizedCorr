######
#' @title corrected
#'
#' @description Calculated the corrected (partial) autocorrelation function
#' @param x a univariate numeric time series.
#' @param lag.max maximum lag at which to calculate the penalized/sample PACF. Defaults to \eqn{10log_{10}(N)} where N is the number of non-missing observations.
#' @param type 'character'. "correlation", "covariance" or "partial".
#' @param na.action function to be called to handle missing values. Default is na.fail, na.pass can be used.
#' @param demean  'logical'. If 'TRUE' (the default), mean(x) is removed prior to estimation.
#' @param lh vector of length 1 (value used for all lags), or length lag.max. Default uses formula in the description.
#' @param ... additional arguments passed to plotting.
#'
#' @return An object of type acf with the following elements:
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
#' # Example for penalized acf estimation
#' corrected(data)
#' }
#' @importFrom stats na.fail
#' @importFrom stats as.ts
#' @importFrom stats toeplitz
#####

corrected=function(x, lag.max = NULL, type = c("correlation", "covariance", 
          "partial"), na.action = na.fail, demean = TRUE, 
          lh=NULL,...){
  x <- na.action(as.ts(x))
  x <- as.matrix(x)
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  sampleT <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  
  if (is.null(lag.max)) 
    lag.max <- floor(10 * (log10(sampleT) - log10(nser)))
  
  h=1:lag.max
  k=max(lag.max,floor(sqrt(sampleT)))
  
  if(is.null(lh)){
    lh=sqrt(log(sampleT)/sampleT)*(k-1+h)/k
  }
  if(length(lh)==1){lh=rep(lh,lag.max)}
  if(length(lh)!=lag.max){stop("lh must either be NULL, length 1, or length lag.max")}
  
  # if pacf then the next line does pacf anyway
  acf=stats::acf(x,lag.max=k,type,plot=F,na.action,demean,...)
  if(type=="partial"){
    tmpacf=acf$acf[1:lag.max,,]
  }# need to take a copy so dimensions are the same as we want the output, original is needed for k lags for bias
  else{
    tmpacf=acf$acf[2:(lag.max+1),,]
  } # need to take a copy so we can remove the lag0 when we do acf and it matches with pacf
  if(nser==1){
    tmpacf=array(tmpacf,dim=c(lag.max,1,1))
  }
  else if(lag.max==1){
    tmpacf=array(tmpacf,dim=c(1,nser,nser))
  }
  
  nserIndexM=matrix(1:nser,ncol=1)
  # bias correction calculation
  j=1:sqrt(acf$n.used)
  if(type=="partial"){
    b=apply(nserIndexM,MARGIN=1,FUN=function(i){
      b=((h+1)*(tmpacf[,i,i]-sign(tmpacf[,i,i])*lh)+(tmpacf[,i,i]+1+h%%2==0))/sampleT
      return(b)
    }) # returns lag.max x nser matrix
  }
  else{
    b=apply(nserIndexM,MARGIN=1,FUN=function(i){
      b=((h+1)*(tmpacf[,i,i]-sign(tmpacf[,i,i])*lh) + 
           (1-tmpacf[,i,i])*(1+2*sum((1-j/sampleT)*acf$acf[j+1,i,i])))/sampleT # +1 as the acf starts at lag0
      return(b)
    }) # returns lag.max x nser matrix
  }
  if(nser==1){
    b=matrix(b,nrow=lag.max,ncol=1)
  }
  else if(lag.max==1){
    b=matrix(b,nrow=1,ncol=nser)
  }
  
  # bias correct if larger than lh, otherwise shrink
  target=apply(nserIndexM,MARGIN=1,FUN=function(i){
    ind=(abs(tmpacf[,i,i])>lh)
    target=(tmpacf[,i,i]+b[,i])*ind+(tmpacf[,i,i]*abs(tmpacf[,i,i])/lh)*!ind
    return(target)
  }) # lag.max x nser
  if(nser==1){
    target=matrix(target,nrow=lag.max,ncol=1)
  }
  else if(lag.max==1){
    target=matrix(target,nrow=1,ncol=nser)
  }
  
  lambda=apply(nserIndexM,MARGIN=1,FUN=function(i){
    ind=(abs(tmpacf[,i,i])>lh)
    a=(lh/tmpacf[,i,i]^2)*!ind 
    b=ind*(sqrt(sampleT)*h*(abs(tmpacf[,i,i])-lh)*(1-lh)/(1-abs(tmpacf[,i,i])))
    lambda=a+b  # stupid R doesn't add the above two if on the same line ?!?
    return(lambda)
  }) # lag.max x nser

  weights=lambda/(1+lambda)
  if(nser==1){
    weights=matrix(weights,nrow=lag.max,ncol=1)
  }
  else if(lag.max==1){
    weights=matrix(weights,nrow=1,ncol=nser)
  }
  
  acfstar=apply(nserIndexM,MARGIN=1,FUN=function(i){
    acfstar=(1-weights[,i])*tmpacf[,i,i] + weights[,i]*target[,i]
    return(acfstar)
  })
  if(nser==1){
    acfstar=matrix(acfstar,nrow=lag.max,ncol=1)
  }
  else if(lag.max==1){
    acfstar=matrix(acfstar,nrow=1,ncol=nser)
  }
  
  # check if NND
  if(type!="partial"){
    acfstar=apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(i){
      tx=toeplitz(c(1,acfstar[,i]))
      ex=min(eigen(tx)$values)
      
      if(ex<=0){
        to=toeplitz(acf$acf[1:(lag.max+1),i,i])
        eo=min(eigen(to)$values)
        
        w=abs(ex)/(eo+abs(ex))
        tx=w*to + (1-w)*tx
        acfstar[,i]=tx[-1,1]
      }
      return(acfstar[,i])
    })
  }
  if(nser==1){
    acfstar=matrix(acfstar,nrow=lag.max,ncol=1)
  }
  else if(lag.max==1){
    acfstar=matrix(acfstar,nrow=1,ncol=nser)
  }
  
  if(type=="partial"){
    acf$acf=acf$acf[1:lag.max,,]
    acf$lag=acf$lag[1:lag.max,,]
    if(nser==1){
      acf$acf=array(acf$acf,dim=c(lag.max,1,1))
      acf$lag=array(acf$lag,dim=c(lag.max,1,1))
    }
    else if(lag.max==1){
      acf$acf=array(acf$acf,dim=c(1,nser,nser))
      acf$lag=array(acf$lag,dim=c(1,nser,nser))
    }
  }
  else{
    acf$acf=acf$acf[1:(lag.max+1),,]
    acf$lag=acf$lag[1:(lag.max+1),,]
    if(nser==1){
      acf$acf=array(acf$acf,dim=c(lag.max+1,1,1))
      acf$lag=array(acf$lag,dim=c(lag.max+1,1,1))
    }
    else if(lag.max==1){
      acf$acf=array(acf$acf,dim=c(2,nser,nser))
      acf$lag=array(acf$lag,dim=c(2,nser,nser))
    }
  }
  
  for(i in 1:nser){
    if(type=="partial"){
      acf$acf[,i,i]=acfstar[,i]
    }
    else{
      acf$acf[-1,i,i]=acfstar[,i]
    }
  }
  acf$lh=lh
  return(acf)
}