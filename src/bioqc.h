/* Common definitions of BioQC */
#ifndef BIOQC_H
#define BIOQC_H

#include <R.h>
#include <Rinternals.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

#endif
