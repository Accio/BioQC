#include "bioqc.h"

/* sortAscending and stat_gini procedures are copied from BIOS library */

static void sortAscending (double x[], int num) {
  /**
     Sorts an array in ascending order
  */
  double tmp;
  int i,k;
  
  for (i=0;i<num-1;i++)
    for (k=i+1;k<num;k++)
      if (x[i] > x[k]) {
        tmp = x[i];
        x[i] = x[k];
        x[k] = tmp;
      }
}

double stat_gini(double x[], int num) {
  /**
     Calculate Gini index.<br>
     Implementation follows R code of package ineq (function Gini(x) )<br>
     Postcondition: x is sorted
     @param[in] x - numbers
     @param[in] num - how many
     @return Gini index
  */ 
  double sum = 0.0;
  double gini = 0.0;
  int i;
  
  sortAscending (x,num);
  for (i=0;i<num;i++) {
    gini += (double)(i+1)*x[i];
    sum += x[i];
  }
  gini = 2.0*(gini/((double)num*sum));
  return gini - 1.0 - (1.0/(double)num);
}

SEXP gini_matrix(SEXP value,
		 SEXP nrowR,
		 SEXP ncolR) {
  double *pmat = REAL(value);
  int nrow = INTEGER(nrowR)[0];
  int ncol = INTEGER(ncolR)[0];
  double rowvec[ncol];
  double curr;
  int i, j, k;
  
  SEXP res;
  PROTECT(res = allocVector(REALSXP,
			    nrow));
  for (i=0; i<nrow; i++) {
    k=0;
    for(j=0; j<ncol; j++) {
      curr=pmat[i+j*nrow];
      if(!ISNA(curr)) {
	rowvec[k++]=curr;
      }
    }
    REAL(res)[i] = stat_gini(rowvec, k);
  }
  UNPROTECT(1);
  return(res);
}
