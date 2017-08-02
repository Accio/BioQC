/* Common definitions of BioQC */
#ifndef BIOQC_H
#define BIOQC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <R.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

// Declarations (prototypes)
extern SEXP gini_numeric(SEXP value, SEXP len);
extern SEXP gini_matrix(SEXP value, SEXP nrowR, SEXP ncolR);
extern SEXP wmw_test(SEXP matrix, SEXP indlist, SEXP rtype);
extern SEXP signed_wmw_test(SEXP matrix, SEXP signedIndList, SEXP rtype);


#ifdef __cplusplus
}
#endif

#endif
