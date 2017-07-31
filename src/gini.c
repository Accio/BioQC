#include <Rinternals.h>
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

double stat_gini_sorted(double xsorted[], int num) {
  /**
     Calculate Gini index.<br>
     Implementation follows R code of package ineq (function Gini(x) )<br>
     Precondition: x is ascendingly sorted and contains no NA
     @param[in] xsorted - numbers ascendingly sorted
     @param[in] num - how many
     @return Gini index
  */ 
  double sum = 0.0;
  double gini = 0.0;
  int i;
  
  for (i=0;i<num;i++) {
    gini += (double)(i+1)*xsorted[i];
    sum += xsorted[i];
  }
  gini = 2.0*(gini/((double)num*sum));
  return gini - 1.0 - (1.0/(double)num);
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
  sortAscending (x,num);
  return(stat_gini_sorted(x, num));
}

extern SEXP gini_numeric(SEXP value, SEXP len) {
  /**
     Calculate Gini index.<br>
     Implementation follows R code of package ineq (function Gini(x) )<br>
     Postcondition: value is sorted ascendingly
     @param[in] value - numbers
     @param[in] len - how many
     @return Gini index
  */ 
  double *ptr = REAL(value);
  int ptrLen = INTEGER(len)[0];
  SEXP res;
  
  PROTECT(res = allocVector(REALSXP, 1));
  REAL(res)[0] = stat_gini_sorted(ptr, ptrLen);
  UNPROTECT(1);
  
  return(res);
}

extern SEXP gini_matrix(SEXP value,
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
