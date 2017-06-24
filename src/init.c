#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP markovchain__list2Mc(SEXP, SEXP, SEXP);
extern SEXP markovchain__matr2Mc(SEXP, SEXP, SEXP, SEXP);
extern SEXP markovchain_canonicForm(SEXP);
extern SEXP markovchain_commclassesKernel(SEXP);
extern SEXP markovchain_commStatesFinder(SEXP);
extern SEXP markovchain_communicatingClasses(SEXP);
extern SEXP markovchain_createSequenceMatrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP markovchain_ctmcFit(SEXP, SEXP, SEXP, SEXP);
extern SEXP markovchain_ExpectedTimeRcpp(SEXP,SEXP);
extern SEXP markovchain_firstpassageKernel(SEXP, SEXP, SEXP);
extern SEXP markovchain_gcd(SEXP, SEXP);
extern SEXP markovchain_generatorToTransitionMatrix(SEXP, SEXP);
extern SEXP markovchain_inferHyperparam(SEXP, SEXP, SEXP);
extern SEXP markovchain_isGen(SEXP);
extern SEXP markovchain_isProb(SEXP);
extern SEXP markovchain_markovchainFit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP markovchain_markovchainListRcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP markovchain_markovchainSequenceParallelRcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP markovchain_markovchainSequenceRcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP markovchain_mcListFitForList(SEXP);
extern SEXP markovchain_multinomCI(SEXP, SEXP, SEXP);
extern SEXP markovchain_multinomialCIForRow(SEXP, SEXP);
extern SEXP markovchain_period(SEXP);
extern SEXP markovchain_predictiveDistribution(SEXP, SEXP, SEXP);
extern SEXP markovchain_priorDistribution(SEXP, SEXP);
extern SEXP markovchain_probabilityatTRCpp(SEXP);
extern SEXP markovchain_recurrentClasses(SEXP);
extern SEXP markovchain_seq2freqProb(SEXP);
extern SEXP markovchain_seq2matHigh(SEXP, SEXP);
extern SEXP markovchain_summaryKernel(SEXP);
extern SEXP markovchain_lexicographicalSort(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"markovchain__list2Mc",                        (DL_FUNC) &markovchain__list2Mc,                         3},
    {"markovchain__matr2Mc",                        (DL_FUNC) &markovchain__matr2Mc,                         4},
    {"markovchain_canonicForm",                     (DL_FUNC) &markovchain_canonicForm,                      1},
    {"markovchain_commclassesKernel",               (DL_FUNC) &markovchain_commclassesKernel,                1},
    {"markovchain_commStatesFinder",                (DL_FUNC) &markovchain_commStatesFinder,                 1},
    {"markovchain_communicatingClasses",            (DL_FUNC) &markovchain_communicatingClasses,             1},
    {"markovchain_createSequenceMatrix",            (DL_FUNC) &markovchain_createSequenceMatrix,             4},
    {"markovchain_ctmcFit",                         (DL_FUNC) &markovchain_ctmcFit,                          4},
    {"markovchain_ExpectedTimeRcpp",                (DL_FUNC) &markovchain_ExpectedTimeRcpp,                 2},
    {"markovchain_firstpassageKernel",              (DL_FUNC) &markovchain_firstpassageKernel,               3},
    {"markovchain_gcd",                             (DL_FUNC) &markovchain_gcd,                              2},
    {"markovchain_generatorToTransitionMatrix",     (DL_FUNC) &markovchain_generatorToTransitionMatrix,      2},
    {"markovchain_inferHyperparam",                 (DL_FUNC) &markovchain_inferHyperparam,                  3},
    {"markovchain_isGen",                           (DL_FUNC) &markovchain_isGen,                            1},
    {"markovchain_isProb",                          (DL_FUNC) &markovchain_isProb,                           1},
    {"markovchain_markovchainFit",                  (DL_FUNC) &markovchain_markovchainFit,                  12},
    {"markovchain_markovchainListRcpp",             (DL_FUNC) &markovchain_markovchainListRcpp,              4},
    {"markovchain_markovchainSequenceParallelRcpp", (DL_FUNC) &markovchain_markovchainSequenceParallelRcpp,  4},
    {"markovchain_markovchainSequenceRcpp",         (DL_FUNC) &markovchain_markovchainSequenceRcpp,          4},
    {"markovchain_mcListFitForList",                (DL_FUNC) &markovchain_mcListFitForList,                 1},
    {"markovchain_multinomCI",                      (DL_FUNC) &markovchain_multinomCI,                       3},
    {"markovchain_multinomialCIForRow",             (DL_FUNC) &markovchain_multinomialCIForRow,              2},
    {"markovchain_period",                          (DL_FUNC) &markovchain_period,                           1},
    {"markovchain_predictiveDistribution",          (DL_FUNC) &markovchain_predictiveDistribution,           3},
    {"markovchain_priorDistribution",               (DL_FUNC) &markovchain_priorDistribution,                2},
    {"markovchain_probabilityatTRCpp",              (DL_FUNC) &markovchain_probabilityatTRCpp,               1},
    {"markovchain_recurrentClasses",                (DL_FUNC) &markovchain_recurrentClasses,                 1},
    {"markovchain_seq2freqProb",                    (DL_FUNC) &markovchain_seq2freqProb,                     1},
    {"markovchain_seq2matHigh",                     (DL_FUNC) &markovchain_seq2matHigh,                      2},
    {"markovchain_summaryKernel",                   (DL_FUNC) &markovchain_summaryKernel,                    1},
    {"markovchain_lexicographicalSort",             (DL_FUNC) &markovchain_lexicographicalSort,              1},
    {NULL, NULL, 0}
};

void R_init_markovchain(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
