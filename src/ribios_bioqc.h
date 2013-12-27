/* Common definitions of ribiosBioQC */
#ifndef RIBIOS_BIOQC_H
#define RIBIOS_BIOQC_H

#include <R.h>
#include <Rinternals.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

#endif
