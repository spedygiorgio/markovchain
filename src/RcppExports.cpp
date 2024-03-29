// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// isGen
bool isGen(NumericMatrix gen);
RcppExport SEXP _markovchain_isGen(SEXP genSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    rcpp_result_gen = Rcpp::wrap(isGen(gen));
    return rcpp_result_gen;
END_RCPP
}
// generatorToTransitionMatrix
NumericMatrix generatorToTransitionMatrix(NumericMatrix gen, bool byrow);
RcppExport SEXP _markovchain_generatorToTransitionMatrix(SEXP genSEXP, SEXP byrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< bool >::type byrow(byrowSEXP);
    rcpp_result_gen = Rcpp::wrap(generatorToTransitionMatrix(gen, byrow));
    return rcpp_result_gen;
END_RCPP
}
// ctmcFit
List ctmcFit(List data, bool byrow, String name, double confidencelevel);
RcppExport SEXP _markovchain_ctmcFit(SEXP dataSEXP, SEXP byrowSEXP, SEXP nameSEXP, SEXP confidencelevelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< bool >::type byrow(byrowSEXP);
    Rcpp::traits::input_parameter< String >::type name(nameSEXP);
    Rcpp::traits::input_parameter< double >::type confidencelevel(confidencelevelSEXP);
    rcpp_result_gen = Rcpp::wrap(ctmcFit(data, byrow, name, confidencelevel));
    return rcpp_result_gen;
END_RCPP
}
// ExpectedTimeRcpp
NumericVector ExpectedTimeRcpp(NumericMatrix x, NumericVector y);
RcppExport SEXP _markovchain_ExpectedTimeRcpp(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(ExpectedTimeRcpp(x, y));
    return rcpp_result_gen;
END_RCPP
}
// probabilityatTRCpp
NumericMatrix probabilityatTRCpp(NumericMatrix y);
RcppExport SEXP _markovchain_probabilityatTRCpp(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(probabilityatTRCpp(y));
    return rcpp_result_gen;
END_RCPP
}
// impreciseProbabilityatTRCpp
NumericVector impreciseProbabilityatTRCpp(S4 C, int i, int t, int s, double error);
RcppExport SEXP _markovchain_impreciseProbabilityatTRCpp(SEXP CSEXP, SEXP iSEXP, SEXP tSEXP, SEXP sSEXP, SEXP errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    rcpp_result_gen = Rcpp::wrap(impreciseProbabilityatTRCpp(C, i, t, s, error));
    return rcpp_result_gen;
END_RCPP
}
// seq2freqProb
NumericVector seq2freqProb(CharacterVector sequence);
RcppExport SEXP _markovchain_seq2freqProb(SEXP sequenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type sequence(sequenceSEXP);
    rcpp_result_gen = Rcpp::wrap(seq2freqProb(sequence));
    return rcpp_result_gen;
END_RCPP
}
// seq2matHigh
NumericMatrix seq2matHigh(CharacterVector sequence, int order);
RcppExport SEXP _markovchain_seq2matHigh(SEXP sequenceSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(seq2matHigh(sequence, order));
    return rcpp_result_gen;
END_RCPP
}
// markovchainSequenceRcpp
CharacterVector markovchainSequenceRcpp(int n, S4 markovchain, CharacterVector t0, bool include_t0);
RcppExport SEXP _markovchain_markovchainSequenceRcpp(SEXP nSEXP, SEXP markovchainSEXP, SEXP t0SEXP, SEXP include_t0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< S4 >::type markovchain(markovchainSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< bool >::type include_t0(include_t0SEXP);
    rcpp_result_gen = Rcpp::wrap(markovchainSequenceRcpp(n, markovchain, t0, include_t0));
    return rcpp_result_gen;
END_RCPP
}
// markovchainListRcpp
List markovchainListRcpp(int n, List object, bool include_t0, CharacterVector t0);
RcppExport SEXP _markovchain_markovchainListRcpp(SEXP nSEXP, SEXP objectSEXP, SEXP include_t0SEXP, SEXP t0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type object(objectSEXP);
    Rcpp::traits::input_parameter< bool >::type include_t0(include_t0SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type t0(t0SEXP);
    rcpp_result_gen = Rcpp::wrap(markovchainListRcpp(n, object, include_t0, t0));
    return rcpp_result_gen;
END_RCPP
}
// markovchainSequenceParallelRcpp
List markovchainSequenceParallelRcpp(S4 listObject, int n, bool include_t0, CharacterVector init_state);
RcppExport SEXP _markovchain_markovchainSequenceParallelRcpp(SEXP listObjectSEXP, SEXP nSEXP, SEXP include_t0SEXP, SEXP init_stateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type listObject(listObjectSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type include_t0(include_t0SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type init_state(init_stateSEXP);
    rcpp_result_gen = Rcpp::wrap(markovchainSequenceParallelRcpp(listObject, n, include_t0, init_state));
    return rcpp_result_gen;
END_RCPP
}
// createSequenceMatrix
NumericMatrix createSequenceMatrix(SEXP stringchar, bool toRowProbs, bool sanitize, CharacterVector possibleStates);
RcppExport SEXP _markovchain_createSequenceMatrix(SEXP stringcharSEXP, SEXP toRowProbsSEXP, SEXP sanitizeSEXP, SEXP possibleStatesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type stringchar(stringcharSEXP);
    Rcpp::traits::input_parameter< bool >::type toRowProbs(toRowProbsSEXP);
    Rcpp::traits::input_parameter< bool >::type sanitize(sanitizeSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type possibleStates(possibleStatesSEXP);
    rcpp_result_gen = Rcpp::wrap(createSequenceMatrix(stringchar, toRowProbs, sanitize, possibleStates));
    return rcpp_result_gen;
END_RCPP
}
// mcListFitForList
List mcListFitForList(List data);
RcppExport SEXP _markovchain_mcListFitForList(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(mcListFitForList(data));
    return rcpp_result_gen;
END_RCPP
}
// _matr2Mc
S4 _matr2Mc(CharacterMatrix matrData, double laplacian, bool sanitize, CharacterVector possibleStates);
RcppExport SEXP _markovchain__matr2Mc(SEXP matrDataSEXP, SEXP laplacianSEXP, SEXP sanitizeSEXP, SEXP possibleStatesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterMatrix >::type matrData(matrDataSEXP);
    Rcpp::traits::input_parameter< double >::type laplacian(laplacianSEXP);
    Rcpp::traits::input_parameter< bool >::type sanitize(sanitizeSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type possibleStates(possibleStatesSEXP);
    rcpp_result_gen = Rcpp::wrap(_matr2Mc(matrData, laplacian, sanitize, possibleStates));
    return rcpp_result_gen;
END_RCPP
}
// _list2Mc
S4 _list2Mc(List data, double laplacian, bool sanitize);
RcppExport SEXP _markovchain__list2Mc(SEXP dataSEXP, SEXP laplacianSEXP, SEXP sanitizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type laplacian(laplacianSEXP);
    Rcpp::traits::input_parameter< bool >::type sanitize(sanitizeSEXP);
    rcpp_result_gen = Rcpp::wrap(_list2Mc(data, laplacian, sanitize));
    return rcpp_result_gen;
END_RCPP
}
// inferHyperparam
List inferHyperparam(NumericMatrix transMatr, NumericVector scale, CharacterVector data);
RcppExport SEXP _markovchain_inferHyperparam(SEXP transMatrSEXP, SEXP scaleSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type transMatr(transMatrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(inferHyperparam(transMatr, scale, data));
    return rcpp_result_gen;
END_RCPP
}
// markovchainFit
List markovchainFit(SEXP data, String method, bool byrow, int nboot, double laplacian, String name, bool parallel, double confidencelevel, bool confint, NumericMatrix hyperparam, bool sanitize, CharacterVector possibleStates);
RcppExport SEXP _markovchain_markovchainFit(SEXP dataSEXP, SEXP methodSEXP, SEXP byrowSEXP, SEXP nbootSEXP, SEXP laplacianSEXP, SEXP nameSEXP, SEXP parallelSEXP, SEXP confidencelevelSEXP, SEXP confintSEXP, SEXP hyperparamSEXP, SEXP sanitizeSEXP, SEXP possibleStatesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    Rcpp::traits::input_parameter< bool >::type byrow(byrowSEXP);
    Rcpp::traits::input_parameter< int >::type nboot(nbootSEXP);
    Rcpp::traits::input_parameter< double >::type laplacian(laplacianSEXP);
    Rcpp::traits::input_parameter< String >::type name(nameSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< double >::type confidencelevel(confidencelevelSEXP);
    Rcpp::traits::input_parameter< bool >::type confint(confintSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type hyperparam(hyperparamSEXP);
    Rcpp::traits::input_parameter< bool >::type sanitize(sanitizeSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type possibleStates(possibleStatesSEXP);
    rcpp_result_gen = Rcpp::wrap(markovchainFit(data, method, byrow, nboot, laplacian, name, parallel, confidencelevel, confint, hyperparam, sanitize, possibleStates));
    return rcpp_result_gen;
END_RCPP
}
// noofVisitsDistRCpp
NumericVector noofVisitsDistRCpp(NumericMatrix matrix, int i, int N);
RcppExport SEXP _markovchain_noofVisitsDistRCpp(SEXP matrixSEXP, SEXP iSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(noofVisitsDistRCpp(matrix, i, N));
    return rcpp_result_gen;
END_RCPP
}
// multinomialCIForRow
NumericMatrix multinomialCIForRow(NumericVector x, double confidencelevel);
RcppExport SEXP _markovchain_multinomialCIForRow(SEXP xSEXP, SEXP confidencelevelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type confidencelevel(confidencelevelSEXP);
    rcpp_result_gen = Rcpp::wrap(multinomialCIForRow(x, confidencelevel));
    return rcpp_result_gen;
END_RCPP
}
// multinomCI
List multinomCI(NumericMatrix transMat, NumericMatrix seqMat, double confidencelevel);
RcppExport SEXP _markovchain_multinomCI(SEXP transMatSEXP, SEXP seqMatSEXP, SEXP confidencelevelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type transMat(transMatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type seqMat(seqMatSEXP);
    Rcpp::traits::input_parameter< double >::type confidencelevel(confidencelevelSEXP);
    rcpp_result_gen = Rcpp::wrap(multinomCI(transMat, seqMat, confidencelevel));
    return rcpp_result_gen;
END_RCPP
}
// commClassesKernel
List commClassesKernel(NumericMatrix P);
RcppExport SEXP _markovchain_commClassesKernel(SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(commClassesKernel(P));
    return rcpp_result_gen;
END_RCPP
}
// communicatingClasses
List communicatingClasses(S4 object);
RcppExport SEXP _markovchain_communicatingClasses(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(communicatingClasses(object));
    return rcpp_result_gen;
END_RCPP
}
// transientStates
CharacterVector transientStates(S4 object);
RcppExport SEXP _markovchain_transientStates(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(transientStates(object));
    return rcpp_result_gen;
END_RCPP
}
// recurrentStates
CharacterVector recurrentStates(S4 object);
RcppExport SEXP _markovchain_recurrentStates(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(recurrentStates(object));
    return rcpp_result_gen;
END_RCPP
}
// recurrentClasses
List recurrentClasses(S4 object);
RcppExport SEXP _markovchain_recurrentClasses(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(recurrentClasses(object));
    return rcpp_result_gen;
END_RCPP
}
// transientClasses
List transientClasses(S4 object);
RcppExport SEXP _markovchain_transientClasses(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(transientClasses(object));
    return rcpp_result_gen;
END_RCPP
}
// reachabilityMatrix
LogicalMatrix reachabilityMatrix(S4 obj);
RcppExport SEXP _markovchain_reachabilityMatrix(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(reachabilityMatrix(obj));
    return rcpp_result_gen;
END_RCPP
}
// isAccessible
bool isAccessible(S4 obj, String from, String to);
RcppExport SEXP _markovchain_isAccessible(SEXP objSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    Rcpp::traits::input_parameter< String >::type from(fromSEXP);
    Rcpp::traits::input_parameter< String >::type to(toSEXP);
    rcpp_result_gen = Rcpp::wrap(isAccessible(obj, from, to));
    return rcpp_result_gen;
END_RCPP
}
// summaryKernel
List summaryKernel(S4 object);
RcppExport SEXP _markovchain_summaryKernel(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(summaryKernel(object));
    return rcpp_result_gen;
END_RCPP
}
// firstpassageKernel
NumericMatrix firstpassageKernel(NumericMatrix P, int i, int n);
RcppExport SEXP _markovchain_firstpassageKernel(SEXP PSEXP, SEXP iSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(firstpassageKernel(P, i, n));
    return rcpp_result_gen;
END_RCPP
}
// firstPassageMultipleRCpp
NumericVector firstPassageMultipleRCpp(NumericMatrix P, int i, NumericVector setno, int n);
RcppExport SEXP _markovchain_firstPassageMultipleRCpp(SEXP PSEXP, SEXP iSEXP, SEXP setnoSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type setno(setnoSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(firstPassageMultipleRCpp(P, i, setno, n));
    return rcpp_result_gen;
END_RCPP
}
// expectedRewardsRCpp
NumericVector expectedRewardsRCpp(NumericMatrix matrix, int n, NumericVector rewards);
RcppExport SEXP _markovchain_expectedRewardsRCpp(SEXP matrixSEXP, SEXP nSEXP, SEXP rewardsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rewards(rewardsSEXP);
    rcpp_result_gen = Rcpp::wrap(expectedRewardsRCpp(matrix, n, rewards));
    return rcpp_result_gen;
END_RCPP
}
// expectedRewardsBeforeHittingARCpp
double expectedRewardsBeforeHittingARCpp(NumericMatrix matrix, int s0, NumericVector rewards, int n);
RcppExport SEXP _markovchain_expectedRewardsBeforeHittingARCpp(SEXP matrixSEXP, SEXP s0SEXP, SEXP rewardsSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< int >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rewards(rewardsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(expectedRewardsBeforeHittingARCpp(matrix, s0, rewards, n));
    return rcpp_result_gen;
END_RCPP
}
// gcd
int gcd(int a, int b);
RcppExport SEXP _markovchain_gcd(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(gcd(a, b));
    return rcpp_result_gen;
END_RCPP
}
// period
int period(S4 object);
RcppExport SEXP _markovchain_period(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(period(object));
    return rcpp_result_gen;
END_RCPP
}
// predictiveDistribution
double predictiveDistribution(CharacterVector stringchar, CharacterVector newData, NumericMatrix hyperparam);
RcppExport SEXP _markovchain_predictiveDistribution(SEXP stringcharSEXP, SEXP newDataSEXP, SEXP hyperparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type stringchar(stringcharSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type newData(newDataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type hyperparam(hyperparamSEXP);
    rcpp_result_gen = Rcpp::wrap(predictiveDistribution(stringchar, newData, hyperparam));
    return rcpp_result_gen;
END_RCPP
}
// priorDistribution
NumericVector priorDistribution(NumericMatrix transMatr, NumericMatrix hyperparam);
RcppExport SEXP _markovchain_priorDistribution(SEXP transMatrSEXP, SEXP hyperparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type transMatr(transMatrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type hyperparam(hyperparamSEXP);
    rcpp_result_gen = Rcpp::wrap(priorDistribution(transMatr, hyperparam));
    return rcpp_result_gen;
END_RCPP
}
// hittingProbabilities
NumericMatrix hittingProbabilities(S4 object);
RcppExport SEXP _markovchain_hittingProbabilities(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type object(objectSEXP);
    rcpp_result_gen = Rcpp::wrap(hittingProbabilities(object));
    return rcpp_result_gen;
END_RCPP
}
// canonicForm
S4 canonicForm(S4 obj);
RcppExport SEXP _markovchain_canonicForm(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(canonicForm(obj));
    return rcpp_result_gen;
END_RCPP
}
// steadyStates
NumericMatrix steadyStates(S4 obj);
RcppExport SEXP _markovchain_steadyStates(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(steadyStates(obj));
    return rcpp_result_gen;
END_RCPP
}
// absorbingStates
CharacterVector absorbingStates(S4 obj);
RcppExport SEXP _markovchain_absorbingStates(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(absorbingStates(obj));
    return rcpp_result_gen;
END_RCPP
}
// isIrreducible
bool isIrreducible(S4 obj);
RcppExport SEXP _markovchain_isIrreducible(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(isIrreducible(obj));
    return rcpp_result_gen;
END_RCPP
}
// isRegular
bool isRegular(S4 obj);
RcppExport SEXP _markovchain_isRegular(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(isRegular(obj));
    return rcpp_result_gen;
END_RCPP
}
// meanAbsorptionTime
NumericVector meanAbsorptionTime(S4 obj);
RcppExport SEXP _markovchain_meanAbsorptionTime(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(meanAbsorptionTime(obj));
    return rcpp_result_gen;
END_RCPP
}
// absorptionProbabilities
NumericMatrix absorptionProbabilities(S4 obj);
RcppExport SEXP _markovchain_absorptionProbabilities(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(absorptionProbabilities(obj));
    return rcpp_result_gen;
END_RCPP
}
// meanFirstPassageTime
NumericMatrix meanFirstPassageTime(S4 obj, CharacterVector destination);
RcppExport SEXP _markovchain_meanFirstPassageTime(SEXP objSEXP, SEXP destinationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type destination(destinationSEXP);
    rcpp_result_gen = Rcpp::wrap(meanFirstPassageTime(obj, destination));
    return rcpp_result_gen;
END_RCPP
}
// meanRecurrenceTime
NumericVector meanRecurrenceTime(S4 obj);
RcppExport SEXP _markovchain_meanRecurrenceTime(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(meanRecurrenceTime(obj));
    return rcpp_result_gen;
END_RCPP
}
// meanNumVisits
NumericMatrix meanNumVisits(S4 obj);
RcppExport SEXP _markovchain_meanNumVisits(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(meanNumVisits(obj));
    return rcpp_result_gen;
END_RCPP
}
// isProb
bool isProb(double prob);
RcppExport SEXP _markovchain_isProb(SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(isProb(prob));
    return rcpp_result_gen;
END_RCPP
}
// isStochasticMatrix
bool isStochasticMatrix(NumericMatrix m, bool byrow);
RcppExport SEXP _markovchain_isStochasticMatrix(SEXP mSEXP, SEXP byrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< bool >::type byrow(byrowSEXP);
    rcpp_result_gen = Rcpp::wrap(isStochasticMatrix(m, byrow));
    return rcpp_result_gen;
END_RCPP
}
// isProbVector
bool isProbVector(NumericVector prob);
RcppExport SEXP _markovchain_isProbVector(SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(isProbVector(prob));
    return rcpp_result_gen;
END_RCPP
}
// checkIsAccesibleMethod
bool checkIsAccesibleMethod(S4 obj);
RcppExport SEXP _markovchain_checkIsAccesibleMethod(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(checkIsAccesibleMethod(obj));
    return rcpp_result_gen;
END_RCPP
}
// approxEqual
bool approxEqual(NumericMatrix a, NumericMatrix b);
RcppExport SEXP _markovchain_approxEqual(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(approxEqual(a, b));
    return rcpp_result_gen;
END_RCPP
}
// isPartition
bool isPartition(List commClasses, CharacterVector states);
RcppExport SEXP _markovchain_isPartition(SEXP commClassesSEXP, SEXP statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type commClasses(commClassesSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type states(statesSEXP);
    rcpp_result_gen = Rcpp::wrap(isPartition(commClasses, states));
    return rcpp_result_gen;
END_RCPP
}
// areHittingProbabilities
bool areHittingProbabilities(NumericMatrix probs, NumericMatrix hitting, bool byrow);
RcppExport SEXP _markovchain_areHittingProbabilities(SEXP probsSEXP, SEXP hittingSEXP, SEXP byrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type hitting(hittingSEXP);
    Rcpp::traits::input_parameter< bool >::type byrow(byrowSEXP);
    rcpp_result_gen = Rcpp::wrap(areHittingProbabilities(probs, hitting, byrow));
    return rcpp_result_gen;
END_RCPP
}
// areMeanNumVisits
bool areMeanNumVisits(NumericMatrix probs, NumericMatrix numVisits, NumericMatrix hitting, bool byrow);
RcppExport SEXP _markovchain_areMeanNumVisits(SEXP probsSEXP, SEXP numVisitsSEXP, SEXP hittingSEXP, SEXP byrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type numVisits(numVisitsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type hitting(hittingSEXP);
    Rcpp::traits::input_parameter< bool >::type byrow(byrowSEXP);
    rcpp_result_gen = Rcpp::wrap(areMeanNumVisits(probs, numVisits, hitting, byrow));
    return rcpp_result_gen;
END_RCPP
}
// recurrentHitting
bool recurrentHitting(List recurrentClasses, NumericMatrix hitting, CharacterVector states, bool byrow);
RcppExport SEXP _markovchain_recurrentHitting(SEXP recurrentClassesSEXP, SEXP hittingSEXP, SEXP statesSEXP, SEXP byrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type recurrentClasses(recurrentClassesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type hitting(hittingSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type states(statesSEXP);
    Rcpp::traits::input_parameter< bool >::type byrow(byrowSEXP);
    rcpp_result_gen = Rcpp::wrap(recurrentHitting(recurrentClasses, hitting, states, byrow));
    return rcpp_result_gen;
END_RCPP
}
// hittingProbsAreOne
bool hittingProbsAreOne(NumericMatrix matrix);
RcppExport SEXP _markovchain_hittingProbsAreOne(SEXP matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(hittingProbsAreOne(matrix));
    return rcpp_result_gen;
END_RCPP
}
// absorbingAreRecurrentClass
bool absorbingAreRecurrentClass(CharacterVector absorbingStates, List recurrentClasses);
RcppExport SEXP _markovchain_absorbingAreRecurrentClass(SEXP absorbingStatesSEXP, SEXP recurrentClassesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type absorbingStates(absorbingStatesSEXP);
    Rcpp::traits::input_parameter< List >::type recurrentClasses(recurrentClassesSEXP);
    rcpp_result_gen = Rcpp::wrap(absorbingAreRecurrentClass(absorbingStates, recurrentClasses));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_markovchain_isGen", (DL_FUNC) &_markovchain_isGen, 1},
    {"_markovchain_generatorToTransitionMatrix", (DL_FUNC) &_markovchain_generatorToTransitionMatrix, 2},
    {"_markovchain_ctmcFit", (DL_FUNC) &_markovchain_ctmcFit, 4},
    {"_markovchain_ExpectedTimeRcpp", (DL_FUNC) &_markovchain_ExpectedTimeRcpp, 2},
    {"_markovchain_probabilityatTRCpp", (DL_FUNC) &_markovchain_probabilityatTRCpp, 1},
    {"_markovchain_impreciseProbabilityatTRCpp", (DL_FUNC) &_markovchain_impreciseProbabilityatTRCpp, 5},
    {"_markovchain_seq2freqProb", (DL_FUNC) &_markovchain_seq2freqProb, 1},
    {"_markovchain_seq2matHigh", (DL_FUNC) &_markovchain_seq2matHigh, 2},
    {"_markovchain_markovchainSequenceRcpp", (DL_FUNC) &_markovchain_markovchainSequenceRcpp, 4},
    {"_markovchain_markovchainListRcpp", (DL_FUNC) &_markovchain_markovchainListRcpp, 4},
    {"_markovchain_markovchainSequenceParallelRcpp", (DL_FUNC) &_markovchain_markovchainSequenceParallelRcpp, 4},
    {"_markovchain_createSequenceMatrix", (DL_FUNC) &_markovchain_createSequenceMatrix, 4},
    {"_markovchain_mcListFitForList", (DL_FUNC) &_markovchain_mcListFitForList, 1},
    {"_markovchain__matr2Mc", (DL_FUNC) &_markovchain__matr2Mc, 4},
    {"_markovchain__list2Mc", (DL_FUNC) &_markovchain__list2Mc, 3},
    {"_markovchain_inferHyperparam", (DL_FUNC) &_markovchain_inferHyperparam, 3},
    {"_markovchain_markovchainFit", (DL_FUNC) &_markovchain_markovchainFit, 12},
    {"_markovchain_noofVisitsDistRCpp", (DL_FUNC) &_markovchain_noofVisitsDistRCpp, 3},
    {"_markovchain_multinomialCIForRow", (DL_FUNC) &_markovchain_multinomialCIForRow, 2},
    {"_markovchain_multinomCI", (DL_FUNC) &_markovchain_multinomCI, 3},
    {"_markovchain_commClassesKernel", (DL_FUNC) &_markovchain_commClassesKernel, 1},
    {"_markovchain_communicatingClasses", (DL_FUNC) &_markovchain_communicatingClasses, 1},
    {"_markovchain_transientStates", (DL_FUNC) &_markovchain_transientStates, 1},
    {"_markovchain_recurrentStates", (DL_FUNC) &_markovchain_recurrentStates, 1},
    {"_markovchain_recurrentClasses", (DL_FUNC) &_markovchain_recurrentClasses, 1},
    {"_markovchain_transientClasses", (DL_FUNC) &_markovchain_transientClasses, 1},
    {"_markovchain_reachabilityMatrix", (DL_FUNC) &_markovchain_reachabilityMatrix, 1},
    {"_markovchain_isAccessible", (DL_FUNC) &_markovchain_isAccessible, 3},
    {"_markovchain_summaryKernel", (DL_FUNC) &_markovchain_summaryKernel, 1},
    {"_markovchain_firstpassageKernel", (DL_FUNC) &_markovchain_firstpassageKernel, 3},
    {"_markovchain_firstPassageMultipleRCpp", (DL_FUNC) &_markovchain_firstPassageMultipleRCpp, 4},
    {"_markovchain_expectedRewardsRCpp", (DL_FUNC) &_markovchain_expectedRewardsRCpp, 3},
    {"_markovchain_expectedRewardsBeforeHittingARCpp", (DL_FUNC) &_markovchain_expectedRewardsBeforeHittingARCpp, 4},
    {"_markovchain_gcd", (DL_FUNC) &_markovchain_gcd, 2},
    {"_markovchain_period", (DL_FUNC) &_markovchain_period, 1},
    {"_markovchain_predictiveDistribution", (DL_FUNC) &_markovchain_predictiveDistribution, 3},
    {"_markovchain_priorDistribution", (DL_FUNC) &_markovchain_priorDistribution, 2},
    {"_markovchain_hittingProbabilities", (DL_FUNC) &_markovchain_hittingProbabilities, 1},
    {"_markovchain_canonicForm", (DL_FUNC) &_markovchain_canonicForm, 1},
    {"_markovchain_steadyStates", (DL_FUNC) &_markovchain_steadyStates, 1},
    {"_markovchain_absorbingStates", (DL_FUNC) &_markovchain_absorbingStates, 1},
    {"_markovchain_isIrreducible", (DL_FUNC) &_markovchain_isIrreducible, 1},
    {"_markovchain_isRegular", (DL_FUNC) &_markovchain_isRegular, 1},
    {"_markovchain_meanAbsorptionTime", (DL_FUNC) &_markovchain_meanAbsorptionTime, 1},
    {"_markovchain_absorptionProbabilities", (DL_FUNC) &_markovchain_absorptionProbabilities, 1},
    {"_markovchain_meanFirstPassageTime", (DL_FUNC) &_markovchain_meanFirstPassageTime, 2},
    {"_markovchain_meanRecurrenceTime", (DL_FUNC) &_markovchain_meanRecurrenceTime, 1},
    {"_markovchain_meanNumVisits", (DL_FUNC) &_markovchain_meanNumVisits, 1},
    {"_markovchain_isProb", (DL_FUNC) &_markovchain_isProb, 1},
    {"_markovchain_isStochasticMatrix", (DL_FUNC) &_markovchain_isStochasticMatrix, 2},
    {"_markovchain_isProbVector", (DL_FUNC) &_markovchain_isProbVector, 1},
    {"_markovchain_checkIsAccesibleMethod", (DL_FUNC) &_markovchain_checkIsAccesibleMethod, 1},
    {"_markovchain_approxEqual", (DL_FUNC) &_markovchain_approxEqual, 2},
    {"_markovchain_isPartition", (DL_FUNC) &_markovchain_isPartition, 2},
    {"_markovchain_areHittingProbabilities", (DL_FUNC) &_markovchain_areHittingProbabilities, 3},
    {"_markovchain_areMeanNumVisits", (DL_FUNC) &_markovchain_areMeanNumVisits, 4},
    {"_markovchain_recurrentHitting", (DL_FUNC) &_markovchain_recurrentHitting, 4},
    {"_markovchain_hittingProbsAreOne", (DL_FUNC) &_markovchain_hittingProbsAreOne, 1},
    {"_markovchain_absorbingAreRecurrentClass", (DL_FUNC) &_markovchain_absorbingAreRecurrentClass, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_markovchain(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
