#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <limits.h>

SEXP spectMat(SEXP x, SEXP n, SEXP window, SEXP offset) {
  if (TYPEOF(x) != REALSXP) error("x must be a double vector");
  if (TYPEOF(window) != REALSXP) error("window must be a double vector");
  if (TYPEOF(offset) != INTSXP) error("offset must be an integer vector");
  if (XLENGTH(n) != 1 || TYPEOF(n) != INTSXP) error("n must be an integer scalar");

  int nRows = INTEGER(n)[0];
  if (nRows == NA_INTEGER || nRows < 0) error("n must be non-negative");

  R_xlen_t nCols = XLENGTH(offset);
  R_xlen_t winSize = XLENGTH(window);
  if (nCols > INT_MAX) error("offset too long");
  if (winSize > INT_MAX) error("window too long");
  if (winSize > nRows) error("length(window) must be <= n");

  R_xlen_t xLen = XLENGTH(x);
  int *ptrOffset = INTEGER(offset);

  for (R_xlen_t i = 0; i < nCols; i++) {
    int off = ptrOffset[i];
    if (off == NA_INTEGER || off < 1)
      error("offset[%lld] must be >= 1", (long long)i + 1);
    R_xlen_t start = (R_xlen_t)off - 1;
    if (start + winSize > xLen)
      error("offset[%lld] + window exceeds length(x)", (long long)i + 1);
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, nRows, (int)nCols));
  double *ptrout = REAL(out);
  double *ptrX = REAL(x);
  double *ptrWindow = REAL(window);

  for (R_xlen_t i = 0; i < nCols; i++) {
    R_xlen_t col = i * (R_xlen_t)nRows;
    R_xlen_t start = (R_xlen_t)ptrOffset[i] - 1;

    for (R_xlen_t j = 0; j < winSize; j++) {
      ptrout[col + j] = ptrX[start + j] * ptrWindow[j];
    }

    if (winSize < (R_xlen_t)nRows) {
      memset(ptrout + col + winSize, 0,
             (size_t)((R_xlen_t)nRows - winSize) * sizeof(double));
    }
  }

  UNPROTECT(1);
  return out;
}
