#include "ribios_bioqc.h"

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

SEXP gini(SEXP value) {
  double g=stat_gini(REAL(value), length(value));

  SEXP res;
  PROTECT(res=allocVector(REALSXP, 1));
  double *pres=REAL(res);
  *pres=g;
  UNPROTECT(1);

  return(res);
}
