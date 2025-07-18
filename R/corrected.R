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
#' @references Gallagher, C., Killick, R., Tan, X. (2024+) Penalized M-estimation 
#' for Autocorrelation. \emph{Submitted.}
#' 
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' # Example for penalized acf estimation
#' corrected(data)
#' }
#' @keywords internal
#' @importFrom stats na.fail
#' @importFrom stats as.ts
#' @importFrom stats toeplitz
#####

corrected=function(x, lag.max = NULL, type = c("correlation", "covariance", 
          "partial"), na.action = na.fail, demean = TRUE, 
          lh=NULL,  
          target=apply(nserIndexM,MARGIN=1,FUN=function(i){
            ind=(abs(tmpacf[,i,i])>lh[,i])
            target=b[,i]*ind
            return(target)
          }),
          lambda=apply(nserIndexM,MARGIN=1,FUN=function(i){
            ind=(abs(tmpacf[,i,i])>lh[,i])
            lambda=(!ind)*h*10*log10(sampleT) *lh[,i]*(lh[,i]-abs(tmpacf[,i,i]))/abs(tmpacf[,i,i])^{3}+ # shrink more aggressively for larger lags
              (ind)*(abs(tmpacf[,i,i])-lh[,i])*(1-lh[,i])/(1-abs(tmpacf[,i,i]))^2*10*log10(sampleT) # movement towards target doesn't depend on lag
            return(lambda)
          }),...){
  x <- na.action(as.ts(x))
  x <- as.matrix(x)
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  sampleT=nrow(x)
  nser=ncol(x)
  
  cov=FALSE
  if(type=="covariance"){
    cov=TRUE
    type=="correlation"
    # set the flag and proceed as if you wanted to calculate correlation
  }

  if(is.null(lag.max)){lag.max <- floor(10 * (log10(sampleT) - log10(nser)))}

  h=1:lag.max
  
  # if pacf then the next line does pacf anyway
  acf=stats::acf(x,lag.max=lag.max,type=type,plot=FALSE,na.action=na.action,demean=demean,...)
  if(type=="partial"){
    tmpacf=acf$acf[1:lag.max,,,drop=FALSE]
  }# need to take a copy so dimensions are the same as we want the output, original is needed for k lags for bias
  else{
    tmpacf=acf$acf[2:(lag.max+1),,,drop=FALSE]
  } # need to take a copy so we can remove the lag0 when we do acf and it matches with pacf
  
  j=1:ceiling(log(sampleT))
  if(type=="partial"){j=j-1}
  
  if(length(lh)==1){lh=matrix(rep(lh,lag.max*nser),nrow=lag.max)} # same value for all series and lags
  else if(length(lh)==lag.max){ # same lh vector for all series
    lh = matrix(lh, nrow = lag.max,ncol=nser)
  }
  else if(length(lh)==nser){ # one value per series, repeat for lags
    lh=matrix(rep(lh,each=lag.max),nrow=lag.max)
  }
  else if(is.null(lh)){ # none supplied so use default
    lh=apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(i){
      el = sqrt(log(sampleT)/(sampleT))
      el=el*sqrt(1-tmpacf[,i,i]^2)
      return(el)
    }) # lh is a matrix lag.max x nser
  }
  else{ # not something we expect so stop
    stop("lh must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
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
    # estimate AR order, if "low" then move toward an inverted PACF rather than bias corrected
    arord=apply(nserIndexM, MARGIN = 1, FUN = function(i) {
      return(ar.penyw(x[,i])$order)
    })
    
    b=apply(nserIndexM,MARGIN=1,FUN=function(i){
      if(arord[i]<(sampleT)^{1/3}){
        # if approximated well by a "low" order AR process
        # using n^{1/3} (had log(n) and log(n)*2 previously), trying to get a balance for larger n as log(n) is too small
        return(acf(x[,i],plot=FALSE,lag.max=lag.max,estimate="invertpacf")$acf[-1,i,i])
      }
      b=tmpacf[,i,i]+(h*tmpacf[,i,i]+(1-tmpacf[,i,i])*(1+2*sum((1-j/sampleT)*tmpacf[j,i,i])))/sampleT
      b = sign(b) * pmin(abs(b), 0.99)
      return(b)
    }) # returns lag.max x nser matrix
  }
  # above doesn't necessarily return a matrix so make one if needed
  if(lag.max==1){b=matrix(b,nrow=1)}
  else if(nser==1){b=matrix(b,ncol=1)}
  
  # bias correct if larger than lh, otherwise shrink
  
  sampleT <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  
  if(!is.numeric(lambda)){stop("'lambda must be numeric")}
  if(!lambda > 0){stop("'lambda' must be positive")}
  if(length(lambda)==1){lambda=matrix(rep(lambda,lag.max*nser),nrow=lag.max)} # same value for all series and lags
  else if(length(lambda)==lag.max){ # same lambda vector for all series
    lambda = matrix(lambda, nrow = lag.max,ncol=nser)
  }
  else if(length(lambda)==nser){ # one value per series, repeat for lags
    lambda=matrix(rep(lambda,each=lag.max),nrow=lag.max)
  }
  else if(is.null(lambda)){ # none supplied so use default
    lambda = apply(nserIndexM,MARGIN=1,FUN=function(i){
      ind=(abs(tmpacf[,i,i])>lh[,i])
      lambda=(!ind)*h*10*log10(sampleT) *lh[,i]*(lh[,i]-abs(tmpacf[,i,i]))/abs(tmpacf[,i,i])^{3}+ # shrink more aggressively for larger lags
        (ind)*(abs(tmpacf[,i,i])-lh[,i])*(1-lh[,i])/(1-abs(tmpacf[,i,i]))^2*10*log10(sampleT) # movement towards target doesn't depend on lag
      return(lambda)
    }) # lambda is a matrix lag.max x nser
  }
  else{ # not something we expect so stop
    stop("lambda must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
  }
  
  
  if(!is.numeric(target)){stop("'target' must be numeric")}
  if(!target < 1){stop("'target' must be between 1 and -1")}
  else if(!target > (-1)){stop("'target' must be between 1 and -1")}
  if(length(target)==1){target=matrix(rep(target,lag.max*nser),nrow=lag.max)} # same value for all series and lags
  else if(length(target)==lag.max){ # same target vector for all series
    target = matrix(target, nrow = lag.max,ncol=nser)
  }
  else if(length(target)==nser){ # one value per series, repeat for lags
    target=matrix(rep(target,each=lag.max),nrow=lag.max)
  }
  else if(is.null(target)){ # none supplied so use default
    target = apply(nserIndexM,MARGIN=1,FUN=function(i){
      ind=(abs(tmpacf[,i,i])>lh[,i])
      target=b[,i]*ind 
      return(target)
    }) # lambda is a matrix lag.max x nser
  }
  else{ # not something we expect so stop
    stop("target must either be NULL, length 1, length lag.max, ncol(x), or a matrix with dimension lag.max x nser.")
  }
  
#  target=apply(nserIndexM,MARGIN=1,FUN=function(i){
#    ind=(abs(tmpacf[,i,i])>lh[,i])
#    target=b[,i]*ind
#    return(target)
#  })  lag.max x nser
#  if(lag.max==1){target=matrix(target,nrow=1)}
#  else if(nser==1){target=matrix(target,ncol=1)}
  
#  lambda=apply(nserIndexM,MARGIN=1,FUN=function(i){
#    ind=(abs(tmpacf[,i,i])>lh[,i])
#    lambda=(!ind)*h*10*log10(sampleT) *lh[,i]*(lh[,i]-abs(tmpacf[,i,i]))/abs(tmpacf[,i,i])^{3}+ # shrink more aggressively for larger lags
#      (ind)*(abs(tmpacf[,i,i])-lh[,i])*(1-lh[,i])/(1-abs(tmpacf[,i,i]))^2*10*log10(sampleT) # movement towards target doesn't depend on lag
#    return(lambda)
#  })  lag.max x nser
  if(lag.max==1){lambda=matrix(lambda,nrow=1)}
  else if(nser==1){lambda=matrix(lambda,ncol=1)}
  
  weights=lambda/(1+lambda) # doesn't drop dimensions

  acfstar=apply(nserIndexM,MARGIN=1,FUN=function(i){
    acfstar=(1-weights[,i])*tmpacf[,i,i] + weights[,i]*target[,i]
            acfstar=acfstar[1:lag.max,drop=FALSE] # needed if lag.max=1
    return(acfstar)
  }) # lag.max x nser

  if(lag.max==1){acfstar=matrix(acfstar,nrow=1)}
  else if(nser==1){acfstar=matrix(acfstar,ncol=1)}
  
  # check if NND
  if(type!="partial"){
    acfstar=apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(i){
        gamma=c(1,acfstar[,i])
        Gamma=toeplitz(gamma)
        ei=eigen(Gamma)
        if(any(ei$values<=0)){ # NND
          if(arord[i]<sampleT^{1/3}){ # modeled by a "low" order AR, move towards that
            # using n^{1/3} (had log(n) and log(n)*2 previously), trying to get a balance for larger n as log(n) is too small
            rr=b[1:lag.max,i]
            R=toeplitz(c(1,rr))
            alpha=min(eigen(R)$values)
            beta=abs(min(ei$values))*(1+1/sampleT)
            cv=beta/(alpha+beta)
            acfstar[,i]=cv*rr+(1-cv)*acfstar[,i]
          }
          else{ # move towards IID 
            R=diag(rep(1,(lag.max+1)))
            alpha=min(eigen(R)$values)
            beta=abs(min(ei$values))*(1+1/sampleT)
            cv=beta/(alpha+beta)
            acfstar[,i]=(1-cv)*acfstar[,i]
          }
          Rfinal=toeplitz(c(1,acfstar[,i]))
        }
        return(acfstar[,i])
    }) # lag.max x nser
  }

  if(type=="partial"){
    acf$acf=acf$acf[1:lag.max,,,drop=FALSE]
    acf$lag=acf$lag[1:lag.max,,,drop=FALSE]
  }
  else{
    acf$acf=acf$acf[1:(lag.max+1),,,drop=FALSE]
    acf$lag=acf$lag[1:(lag.max+1),,,drop=FALSE]
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
