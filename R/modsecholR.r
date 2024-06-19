secholC=function(B){
  # R wrapper for C function contained in modsechol.c and called modchl
  
  #   int *ndim;				/* dimension of A */
  #   int *n;					/* largest dimension of matrix that will be used */
  #   double *A;				/* n*n symmetric matrix (only lower triangular portion of A, including the main diagonal, is used) */
  #   double *mcheps;	/* machine precision */
  #   double *tau1;		/* tolerance used for determining when to switch phase 2 */
  #   double *tau2;		/* tolerance used for determining the maximum condition number of the final 2X2 submatrix. */
  
#  dyn.load('modsechol.so')
	mcheps=.Machine$double.eps
	tau1 = mcheps* (1/3)
	tau2 = mcheps* (1/3)
  
  storage.mode(B)='double'
  answer=.C('modchl',ndim=as.integer(nrow(B)),n=as.integer(nrow(B)),A=B,mcheps=as.numeric(mcheps),
            tau1=as.numeric(tau1),tau2=as.numeric(tau2),PACKAGE="PenalizedCorr")
  return(answer$A)
}

