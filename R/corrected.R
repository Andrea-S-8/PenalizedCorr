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
  j=1:floor(sqrt(sampleT))
  if(type=="partial"){j=j-1}
  if(is.null(lh)){
    lh=apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(i){
      el=sqrt(log(sampleT)/(sampleT)*(1-tmpacf[,i,i]^2))*((k-1+h)/k)
      if(lag.max>1){el=(var((acf$acf[j+1,i,i]))/mean(acf$acf[j+1,i,i]^2))*el}
      return(el)
    }) # lh is a matrix lag.max x nser
  }
  if(length(lh)==1){lh=matrix(rep(lh,lag.max*nser),nrow=lag.max)}
  if(length(lh)==lag.max){lh=matrix(rep(lh,nser),nrow=lag.max)}
  if(any(dim(lh)!=c(lag.max,nser))){
    stop("lh must either be NULL, length 1, length lag.max or a matrix with dimension lag.max x nser.")
  }
  
  nserIndexM=matrix(1:nser,ncol=1)
  # bias correction calculation
  if(type=="partial"){
    b=apply(nserIndexM,MARGIN=1,FUN=function(i){
      b=tmpacf[,i,i]+((h+1)*tmpacf[,i,i]+(tmpacf[,i,i]+1+h%%2==0))/sampleT
              b=sign(b)*pmin(abs(b),0.99)
      return(b)
    }) # returns lag.max x nser matrix
  }
  else{
    b=apply(nserIndexM,MARGIN=1,FUN=function(i){
      b=tmpacf[,i,i]+((h+1)*tmpacf[,i,i] + 
            (1-tmpacf[,i,i])*(1+2*sum((1-j/sampleT)*acf$acf[j+1,i,i])))/sampleT # +1 as the acf starts at lag0
              b=sign(b)*pmin(abs(b),0.99)
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
    ind=(abs(b[,i])>lh[,i])
    target=b[,i]*ind
    return(target)
  }) # lag.max x nser
  if(nser==1){
    target=matrix(target,nrow=lag.max,ncol=1)
  }
  else if(lag.max==1){
    target=matrix(target,nrow=1,ncol=nser)
  }
  
  lambda=apply(nserIndexM,MARGIN=1,FUN=function(i){
    ind=(abs(b[,i])>lh[,i])
    lambda=(!ind)*h*lh[,i]*(lh[,i]-abs(b[,i]))/abs(b[,i])^{3}+(ind)*h*(abs(b[,i])-lh[,i])*(1-lh[,i])/(1-abs(b[,i]))^2*sqrt(sampleT)
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
        gamma=c(1,acfstar[,i])
        Gamma <- toeplitz(gamma)
        # Compute eigenvalue decomposition
        ei <- eigen(Gamma)
        if(any(ei$values<=0)){ # original is NND
          # Shrink eigenvalues
          d <- pmax(ei$values, 10 / sampleT)
          # Construct new covariance matrix
          Gamma2 <- ei$vectors %*% diag(d) %*% t(ei$vectors)
          Gamma2 <- Gamma2 / mean(d)
          # Estimate new ACF
          d <- row(Gamma2) - col(Gamma2)
          for (j in 2:(lag.max+1))
            gamma[j] <- mean(Gamma2[d == (j - 1)]) 
          
          acfstar[,i]=gamma[-1]
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
