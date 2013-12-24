#include "ribios_bioqc.h"

SEXP gini(SEXP value) {
  double g=stat_gini(REAL(value), length(value));

  SEXP res;
  PROTECT(res=allocVector(REALSXP, 1));
  double *pres=REAL(res);
  *pres=g;
  UNPROTECT(1);

  return(res);
}
