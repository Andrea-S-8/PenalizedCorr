#include <R.h> 
#include <Rmath.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

void modchl(int *ndim,int *n,double *A,double *mcheps,double *tau1,double *tau2){
/*	int *ndim;				 dimension of A */
/*	int *n;					 largest dimension of matrix that will be used */
/*	double *A;				 n*n symmetric matrix (only lower triangular portion of A, including the main diagonal, is used) */
/*	double *mcheps;	 machine precision */
/*	double *tau1;		 tolerance used for determining when to switch phase 2 */
/*	double *tau2;		 tolerance used for determining the maximum condition number of the final 2X2 submatrix. */
  
		/* performs a modified cholesky factorization of the form (Ptranspose)AP  + E = L(Ltranspose),
       where L is stored in the lower triangle of the original matrix A. */

		/* create all the variables we are going to need */
		register int j;        /*- current iteration number */
		register int iming;    /*- index of the row with the min. of the neg. lower Gersch. bounds */
		register int imaxd;    /*- index of the row with the maximum diag. element */
		register int i,itemp,jp1,k;  /*- temporary integer variables */
		double delta=0.0;    /*- amount to add to Ajj at the jth iteration */
		double gamma=0.0;    /*- the maximum diagonal element of the original matrix A. */
		double normj;    /*- the 1 norm of A(colj), rows j+1 --> n. */
		double ming=0.0;     /*- the minimum of the neg. lower Gersch. bounds */
		double maxd;     /*- the maximum diagonal element */
		double taugam;		/*- tau1 * gamma */
		int phase1=1;      /*- logical, true if in phase1, otherwise false */
		double delta1,temp,jdmin,tdmin,tempjj;  /*- temporary double precision vars. */
		double g[*n];			/* - n*1 work array */
		int P[*n];			/* - a record of how the rows and columns of the matrix were permuted while performing the decomposition */
		double E[*n];					/* - n*1 array, the ith element is the amount added to the diagonal of A at the ith iteration */

//		######### Initialise ###########
		for(i=1;i<(*n+1);++i){
			g[i-1]=0;
			P[i-1]=i;
			E[i-1]=0;
		}
//     Find the maximum magnitude of the diagonal elements.
//     If any diagonal element is negative, then phase1 is false.

		for(i=1;i<(*n-1);++i){
			if(gamma<fabs(*(A+i-1+(i-1)*(*n)))){ gamma=fabs(*(A+i-1+(i-1)*(*n))); }
			if (*(A+i-1+(i-1)*(*n)) < 0.0){ phase1=0; }
		}

		taugam = *tau1 * gamma;

//     If not in phase1, then calculate the initial gerschgorin bounds
//     needed for the start of phase2.

		if(phase1==0){
			j=1;
			for(i=j;i<(*n+1);++i){
				double offrow = 0.0;
				if(i!=j){
					for(k=j;k<i;++k){
// +1 on the memory address in C means add one to the column in R (works right to left in dimensions)
						offrow += fabs(*(A+(i-1)+(k-1)*(*n)));
					}
				}
				if(i!=*n){
					for(k=(i+1);k<(*n+1);++k){
						offrow += fabs(*(A+(i-1)*(*n)+(k-1)));
					}
				}
				g[i-1] = offrow - *(A+i-1+(i-1)*(*n));
			}
		}

//     check for n=1
		if(*n==1){
			delta = (*tau2 * abs(*A)) - *A ;
			if (delta > 0){ E[1] = delta; }
			if (*A == 0){ E[1] = *tau2; }
			*A=sqrt(*A+E[1]);
		}



		for(j=1;j<(*n);++j){
//        PHASE 1
			if (phase1==1){
//           Find index of maximum diagonal element A(i,i) where i>=j
				maxd = *(A+(j-1)*(*n)+(j-1));
				imaxd = j;
				for(i=(j+1);i<(*n+1);++i){
					if (maxd < *(A+(i-1)*(*n)+(i-1))){
						maxd = *(A+(i-1)*(*n)+(i-1));
						imaxd = i;
					}
				}

//           Pivot to the top the row and column with the max diag
				if (imaxd != j){
//              Swap row j with row of max diag
					for(i=1;i<j;++i){
						temp = *(A+j-1+(i-1)*(*n));
						*(A+j-1+(i-1)*(*n)) = *(A+imaxd-1+(i-1)*(*n));
						*(A+imaxd-1+(i-1)*(*n)) = temp;
					}
//              Swap colj and row maxdiag between j and maxdiag
					for(i=(j+1);i<imaxd;++i){
						temp = *(A+i-1+(j-1)*(*n));
						*(A+i-1+(j-1)*(*n)) = *(A+imaxd-1+(i-1)*(*n));
						*(A+imaxd-1+(i-1)*(*n)) = temp;
					}
//              Swap column j with column of max diag
					if(imaxd!=*n){
						for(i=(imaxd+1);i<(*n+1);++i){
							temp = *(A+i-1+(j-1)*(*n));
							*(A+i-1+(j-1)*(*n)) = *(A+i-1+(imaxd-1)*(*n));
							*(A+i-1+(imaxd-1)*(*n)) = temp;
						}
					}
//              Swap diag elements
					temp = *(A+j-1+(j-1)*(*n));
					*(A+j-1+(j-1)*(*n)) = *(A+imaxd-1+(imaxd-1)*(*n));
					*(A+imaxd-1+(imaxd-1)*(*n)) = temp;
//              Swap elements of the permutation vector
					itemp = P[j];
					P[j] = P[imaxd];
					P[imaxd] = itemp;
				}

//           Check to see whether the normal cholesky update for this
//           iteration would result in a positive diagonal,
//           and if not then switch to phase 2.
				jp1 = j+1;
				tempjj=*(A+j-1+(j-1)*(*n));
				if (tempjj>0){
					jdmin=*(A+jp1-1+(jp1-1)*(*n));
					for(i=jp1;i<(*n+1);++i){
						temp = *(A+i-1+(j-1)*(*n)) * (*(A+i-1+(j-1)*(*n)) / tempjj);
						tdmin = *(A+i-1+(i-1)*(*n)) - temp;
						if(jdmin>tdmin){jdmin = tdmin;}
					}
					if (jdmin < taugam){ phase1 = 0;}
				}
				else{
					phase1 = 0;
				}

				if (phase1==1){
//              Do the normal cholesky update if still in phase 1
					*(A+j-1+(j-1)*(*n)) = sqrt(*(A+j-1+(j-1)*(*n)));
					tempjj = *(A+j-1+(j-1)*(*n));
					for(i=jp1;i<(*n+1);++i){
						*(A+i-1+(j-1)*(*n)) = *(A+i-1+(j-1)*(*n)) / tempjj;
					}
					for(i=jp1;i<(*n+1);++i){
						temp=*(A+i-1+(j-1)*(*n));
						for(k=jp1;k<(i+1);++k){
							*(A+i-1+(k-1)*(*n)) -= (temp * (*(A+k-1+(j-1)*(*n))));
						}
					}
					if (j == (*n-1)){ *(A+*n-1+(*n-1)*(*n))=sqrt(*(A+*n-1+(*n-1)*(*n))); }
				}
				else{
//              Calculate the negatives of the lower gerschgorin bounds
					int ig,kg;
					double offrow;
					for(ig=j;ig<(*n+1);++ig){
						offrow=0.0;
						if(ig!=j){
							for(kg=j;kg<ig;++kg){
								offrow += abs(*(A+ig-1+(kg-1)*(*n)));
							}
						}
						if(ig!=*n){
							for(kg=(ig+1);kg<(*n+1);++kg){
								offrow += abs(*(A+kg-1+(ig-1)*(*n)));
							}
						}
            g[ig] = offrow - *(A+ig-1+(ig-1)*(*n));
					}

				}
			}

//        PHASE 2

			else if (phase1==0){
				if (j !=(*n-1)){
//              Find the minimum negative gershgorin bound
					iming=j;
					ming = g[j];
					for(i=(j+1);i<(*n+1);++i){
						if (ming > g[i]){
							ming = g[i];
							iming = i;
						}
					}
//               Pivot to the top the row and column with the
//              minimum negative gerschgorin bound
					if (iming != j){
//                  Swap row j with row of min gersch bound
						for(i=1;i<j;++i){
							temp = *(A+j-1+(i-1)*(*n));
							*(A+j-1+(i-1)*(*n)) = *(A+iming-1+(i-1)*(*n));
							*(A+iming-1+(i-1)*(*n)) = temp;
						}
//                  Swap colj with row iming from j to iming
						for(i=(j+1);i<iming;++i){
							temp = *(A+i-1+(j-1)*(*n));
							*(A+i-1+(j-1)*(*n)) = *(A+iming-1+(i-1)*(*n));
							*(A+iming-1+(i-1)*(*n)) = temp;
						}
//                 Swap column j with column of min gersch bound
						if(iming!=*n){
							for(i=(iming+1);i<(*n+1);++i){
								temp = *(A+i-1+(j-1)*(*n));
								*(A+i-1+(j-1)*(*n)) = *(A+i-1+(iming-1)*(*n));
								*(A+i-1+(iming-1)*(*n)) = temp;
							}
						}
//                 Swap diagonal elements
						temp = *(A+j-1+(j-1)*(*n));
						*(A+j-1+(j-1)*(*n)) = *(A+iming-1+(iming-1)*(*n));
						*(A+iming-1+(iming-1)*(*n)) = temp;
//                 Swap elements of the permutation vector
						itemp = P[j];
						P[j] = P[iming];
						P[iming] = itemp;
//                 Swap elements of the negative gerschgorin bounds vecto
						temp = g[j];
						g[j] = g[iming];
						g[iming] = temp;
					}

//              Calculate delta and add to the diagonal.
//              delta=max{0,-A(j,j) + max{normj,taugam},delta_previous}
//              where normj=sum of |A(i,j)|,for i=1,n,
//              delta_previous is the delta computed at the previous iter
//              and taugam is tau1*gamma.

					normj = 0.0;
					for(i=(j+1);i<(*n+1);++i){
						normj += abs(*(A+i-1+(j-1)*(*n)));
					}
					temp = normj;
					if(temp<taugam){temp=taugam;}
					delta1 = temp - *(A+j-1+(j-1)*(*n));
					if(delta1<0){delta1=0;}
					if(delta<delta1){delta=delta1;}
					E[j] =  delta;
					*(A+j-1+(j-1)*(*n)) += E[j];

//              Update the gerschgorin bound estimates
//              (note: g[i] is the negative of the
//               Gerschgorin lower bound.)
					if (*(A+j-1+(j-1)*(*n)) != normj){
						temp = (normj/ *(A+j-1+(j-1)*(*n))) - 1.0;
						for(i=(j+1);i<(*n+1);++i){
							g[i] += abs(*(A+i-1+(j-1)*(*n))) * temp;
						}
					}
//              Do the cholesky update
					*(A+j-1+(j-1)*(*n)) = sqrt(*(A+j-1+(j-1)*(*n)));
					tempjj = *(A+j-1+(j-1)*(*n));
					for(i=(j+1);i<(*n+1);++i){
						*(A+i-1+(j-1)*(*n)) = *(A+i-1+(j-1)*(*n)) / tempjj;
					}
					for(i=(j+1);i<(*n+1);++i){
						temp = *(A+i-1+(j-1)*(*n));
						for(k=(j+1);k<(i+1);++k){
							*(A+i-1+(k-1)*(*n)) -= (temp * *(A+k-1+(j-1)*(*n)));
						}
					}
				}
				else{
					double t1,t2,t3;
					double lmbd1,lmbd2,lmbdhi,lmbdlo;
//     Find eigenvalues of final 2 by 2 submatrix
					t1 = *(A+*n-2+(*n-2)*(*n)) + *(A+*n-1+(*n-1)*(*n));
					t2 = *(A+*n-2+(*n-2)*(*n)) - *(A+*n-1+(*n-1)*(*n));
					t3 = sqrt(t2*t2 + 4.0*(*(A+*n-1+(*n-2)*(*n)))*(*(A+*n-1+(*n-2)*(*n))));
					lmbd1 = (t1 - t3)/2;
					lmbd2 = (t1 + t3)/2;
					lmbdhi=lmbd1;
					if(lmbdhi<lmbd2){lmbdhi=lmbd2;}
					lmbdlo=lmbd1;
					if(lmbdlo>lmbd2){lmbdlo=lmbd2;}

//     Find delta such that:
//     1.  the l2 condition number of the final
//     2X2 submatrix + delta*I <= tau2
//     2. delta >= previous delta,
//     3. lmbdlo + delta >= tau2 * gamma,
//     where lmbdlo is the smallest eigenvalue of the final
//     2X2 submatrix

					delta1=(lmbdhi-lmbdlo)/(1.0-*tau2);
					if(delta1<gamma){delta1= gamma;}
					delta1= *tau2 * delta1 - lmbdlo;
					if(delta<0){delta = 0;}
					if(delta<delta1){delta=delta1;}

					if (delta > 0){
						*(A+*n-2+(*n-2)*(*n)) += delta;
						*(A+*n-1+(*n-1)*(*n)) += delta;
						E[*n-1] = delta;
						E[*n] = delta;
					}

//     Final update
					*(A+*n-2+(*n-2)*(*n)) = sqrt(*(A+*n-2+(*n-2)*(*n)));
					*(A+*n-1+(*n-2)*(*n)) = *(A+*n-1+(*n-2)*(*n))/ *(A+*n-2+(*n-2)*(*n));
					*(A+*n-1+(*n-1)*(*n)) -= (*(A+*n-1+(*n-2)*(*n)))* (*(A+*n-1+(*n-2)*(*n)));
					*(A+*n-1+(*n-1)*(*n)) = sqrt(*(A+*n-1+(*n-1)*(*n)));

				}
			}
		}
}

