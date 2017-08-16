#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _SuperRanker_overlap(SEXP);
extern SEXP _SuperRanker_sracppfull(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_SuperRanker_overlap",    (DL_FUNC) &_SuperRanker_overlap,    1},
    {"_SuperRanker_sracppfull", (DL_FUNC) &_SuperRanker_sracppfull, 2},
    {NULL, NULL, 0}
};

void R_init_SuperRanker(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
