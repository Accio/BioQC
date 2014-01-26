#include <R_ext/Rdynload.h>
#include "bioqc.h"

static const R_CallMethodDef callMethods[] = {
  CALLMETHOD_DEF("gini", 1),
  CALLMETHOD_DEF("read_gmt", 1),
  {NULL, NULL, 0}
};

void R_init_ribiosBioQC(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
