invertpacf=function(x, lag.max = NULL, type = c("correlation", "covariance", "partial"), 
                   na.action = na.fail, demean = TRUE, penalized=TRUE,lh=NULL,
                   estarorder=TRUE,...){
  x <- na.action(as.ts(x))
  x <- as.matrix(x)
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  sampleT=nrow(x)
  nser=ncol(x)
  if (is.null(lag.max)) 
    lag.max <- floor(10 * (log10(sampleT) - log10(nser)))
  if(is.matrix(lh)){
    if(dim(lh)[1]!=lag.max){
      stop("lh is a matrix with the incorrect dimension, must have nrow(lh)=lag.max.")
    }
    else if(dim(lh)[2]!=nser){
      stop("lh is a matrix with the incorrect dimension, must have ncol(lh)=nser.")
    }
  }
  
  # First get the approximate AR order by AIC
  xpacf=pacf(x,lag.max=lag.max,plot=F,na.action=na.action,penalized=penalized,lh=lh)
  
  xarorder=lag.max # just the maximum lag if estarorder=FALSE
  if(estarorder){
    AICpen <- apply(matrix(1:lag.max,ncol=1), MARGIN=1,FUN=function(i){
      sampleT*log(apply(x,MARGIN=2,FUN=var)* # nser length of variances
                    apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(j){
                      prod(1-xpacf$acf[1:i,j,j]^2) # nser length of products
                    })) + 2*i
    })
    if(nser==1){AICpen=matrix(AICpen,nrow=1)}
    AICpen <- cbind(sampleT * log(apply(x,MARGIN=2,FUN=var)),AICpen)
    
    xarorder <- apply(AICpen,MARGIN=1,FUN=which.min)-1
    for(i in 1:nser){ # zero the pacf above the fitted ar lag
      if(xarorder[i]<lag.max){
        xpacf$acf[(xarorder[i]+1):lag.max,i,i]=0
      }
    }
  }
  
  
  # Then calculate the AR coefficients for that order
  arcoef=apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(i){
    if(xarorder[i]==0){
      return(NULL)
    }
    if(is.matrix(lh)){
      arcoef=DLpencoef(x[,i], lag.max = xarorder[i], na.action=na.action,penalized=penalized,lh=lh[,i],return.mat=TRUE)
    }
    else{
      arcoef=DLpencoef(x[,i], lag.max = xarorder[i], na.action=na.action,penalized=penalized,lh=lh,return.mat=TRUE)
    }
    if(xarorder[i]==1){
      arcoef=matrix(arcoef,nrow=1)
    }
    return(arcoef)
  },simplify=FALSE) # list of length nser with each element being xarorder[i] x xarorder[i]
  
  # Then convert those into the acf
  acf=stats::acf(x,lag.max=lag.max,type=type,plot=FALSE,na.action=na.action,demean=demean)
  xacf=apply(matrix(1:nser,ncol=1),MARGIN=1,FUN=function(i){
    if(xarorder[i]==0){ # Independence so no AR terms
      xacf=rep(0,lag.max)
      return(xacf)
    }
    xacf=xpacf$acf[,i,i,drop=FALSE] #initialize
    prodacf=cumprod((1-xpacf$acf[1:(lag.max-1),i,i]^2)) * xpacf$acf[2:lag.max,i,i]
    if(xarorder[i]!=1){
      for(j in 2:min(lag.max,xarorder[i])){
        xacf[j]=prodacf[j-1] + arcoef[[i]][1:(j-1),j-1]%*%xacf[(j-1):1]
      }
    }
    if((xarorder[i]>=lag.max)){return(xacf)}
    for(j in (xarorder[i]+1):lag.max){
      xacf[j]=arcoef[[i]][,xarorder[i]]%*%xacf[(j-1):(j-xarorder[i])]
    }
    return(xacf)
  })
  
  for(i in 1:nser){
    acf$acf[-1,i,i]=xacf[,i]
  }
  return(acf)
}