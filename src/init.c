#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern SEXP spectMat(SEXP x, SEXP n, SEXP window, SEXP offset);

static const R_CallMethodDef CallEntries[] = {
  {"spectMat", (DL_FUNC) &spectMat, 4},
  {NULL, NULL, 0}
};

void R_init_Ruido(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
