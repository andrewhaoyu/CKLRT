// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// Eigen_C
List Eigen_C(NumericMatrix As);
RcppExport SEXP _CKLRT_Eigen_C(SEXP AsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type As(AsSEXP);
    rcpp_result_gen = Rcpp::wrap(Eigen_C(As));
    return rcpp_result_gen;
END_RCPP
}
// Eigen_C_value
Eigen::VectorXd Eigen_C_value(NumericMatrix As);
RcppExport SEXP _CKLRT_Eigen_C_value(SEXP AsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type As(AsSEXP);
    rcpp_result_gen = Rcpp::wrap(Eigen_C_value(As));
    return rcpp_result_gen;
END_RCPP
}
// MatMult_C
SEXP MatMult_C(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _CKLRT_MatMult_C(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(MatMult_C(A, B));
    return rcpp_result_gen;
END_RCPP
}
// Sum_C
double Sum_C(NumericVector AA);
RcppExport SEXP _CKLRT_Sum_C(SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(Sum_C(AA));
    return rcpp_result_gen;
END_RCPP
}
// ColSum_C
NumericVector ColSum_C(NumericMatrix AA);
RcppExport SEXP _CKLRT_ColSum_C(SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(ColSum_C(AA));
    return rcpp_result_gen;
END_RCPP
}
// MatrixRowMax_C
NumericVector MatrixRowMax_C(NumericMatrix AA);
RcppExport SEXP _CKLRT_MatrixRowMax_C(SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(MatrixRowMax_C(AA));
    return rcpp_result_gen;
END_RCPP
}
// Elementwisesquare_C
NumericVector Elementwisesquare_C(NumericMatrix AA);
RcppExport SEXP _CKLRT_Elementwisesquare_C(SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(Elementwisesquare_C(AA));
    return rcpp_result_gen;
END_RCPP
}
// VecMultMat_C
NumericVector VecMultMat_C(Eigen::VectorXd A, Eigen::MatrixXd B);
RcppExport SEXP _CKLRT_VecMultMat_C(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(VecMultMat_C(A, B));
    return rcpp_result_gen;
END_RCPP
}
// Vecplus_C
NumericVector Vecplus_C(Eigen::VectorXd A, Eigen::VectorXd B);
RcppExport SEXP _CKLRT_Vecplus_C(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(Vecplus_C(A, B));
    return rcpp_result_gen;
END_RCPP
}
// ColSumtwomatrix_C
NumericVector ColSumtwomatrix_C(NumericMatrix AA, NumericMatrix BB);
RcppExport SEXP _CKLRT_ColSumtwomatrix_C(SEXP AASEXP, SEXP BBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type AA(AASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type BB(BBSEXP);
    rcpp_result_gen = Rcpp::wrap(ColSumtwomatrix_C(AA, BB));
    return rcpp_result_gen;
END_RCPP
}
// ifelsetest_C
NumericVector ifelsetest_C(NumericVector x);
RcppExport SEXP _CKLRT_ifelsetest_C(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(ifelsetest_C(x));
    return rcpp_result_gen;
END_RCPP
}
// MatrixPlus_C
SEXP MatrixPlus_C(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _CKLRT_MatrixPlus_C(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(MatrixPlus_C(A, B));
    return rcpp_result_gen;
END_RCPP
}
// NumxMatrix_C
SEXP NumxMatrix_C(double A, Eigen::MatrixXd B);
RcppExport SEXP _CKLRT_NumxMatrix_C(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(NumxMatrix_C(A, B));
    return rcpp_result_gen;
END_RCPP
}
// LR0_fixRho_C
NumericMatrix LR0_fixRho_C(NumericVector LamdasR, NumericVector muR, NumericMatrix w1R, NumericMatrix w2R, int nminuspx);
RcppExport SEXP _CKLRT_LR0_fixRho_C(SEXP LamdasRSEXP, SEXP muRSEXP, SEXP w1RSEXP, SEXP w2RSEXP, SEXP nminuspxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type LamdasR(LamdasRSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type muR(muRSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w1R(w1RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w2R(w2RSEXP);
    Rcpp::traits::input_parameter< int >::type nminuspx(nminuspxSEXP);
    rcpp_result_gen = Rcpp::wrap(LR0_fixRho_C(LamdasR, muR, w1R, w2R, nminuspx));
    return rcpp_result_gen;
END_RCPP
}
// doubleloop
NumericMatrix doubleloop(NumericMatrix K1R, NumericMatrix K2R, NumericMatrix P0R, NumericMatrix AR, NumericMatrix U1R, NumericMatrix wR, NumericVector LamdasR, int nminuspx, NumericVector all_rho, NumericMatrix LR0_allRhoR);
RcppExport SEXP _CKLRT_doubleloop(SEXP K1RSEXP, SEXP K2RSEXP, SEXP P0RSEXP, SEXP ARSEXP, SEXP U1RSEXP, SEXP wRSEXP, SEXP LamdasRSEXP, SEXP nminuspxSEXP, SEXP all_rhoSEXP, SEXP LR0_allRhoRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type K1R(K1RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K2R(K2RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P0R(P0RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type AR(ARSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type U1R(U1RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wR(wRSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LamdasR(LamdasRSEXP);
    Rcpp::traits::input_parameter< int >::type nminuspx(nminuspxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type all_rho(all_rhoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type LR0_allRhoR(LR0_allRhoRSEXP);
    rcpp_result_gen = Rcpp::wrap(doubleloop(K1R, K2R, P0R, AR, U1R, wR, LamdasR, nminuspx, all_rho, LR0_allRhoR));
    return rcpp_result_gen;
END_RCPP
}
// LR0_fixRho_LRT_C
NumericMatrix LR0_fixRho_LRT_C(NumericVector LamdasR, NumericVector muR, NumericMatrix w1R, NumericMatrix w2R, int nminuspx, NumericVector xiR);
RcppExport SEXP _CKLRT_LR0_fixRho_LRT_C(SEXP LamdasRSEXP, SEXP muRSEXP, SEXP w1RSEXP, SEXP w2RSEXP, SEXP nminuspxSEXP, SEXP xiRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type LamdasR(LamdasRSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type muR(muRSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w1R(w1RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w2R(w2RSEXP);
    Rcpp::traits::input_parameter< int >::type nminuspx(nminuspxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xiR(xiRSEXP);
    rcpp_result_gen = Rcpp::wrap(LR0_fixRho_LRT_C(LamdasR, muR, w1R, w2R, nminuspx, xiR));
    return rcpp_result_gen;
END_RCPP
}
// doubleloop_LRT
NumericMatrix doubleloop_LRT(NumericMatrix K1R, NumericMatrix K2R, NumericMatrix P0R, NumericMatrix AR, NumericMatrix U1R, NumericMatrix wR, NumericVector LamdasR, int nminuspx, NumericVector all_rho, NumericMatrix LR0_allRhoR);
RcppExport SEXP _CKLRT_doubleloop_LRT(SEXP K1RSEXP, SEXP K2RSEXP, SEXP P0RSEXP, SEXP ARSEXP, SEXP U1RSEXP, SEXP wRSEXP, SEXP LamdasRSEXP, SEXP nminuspxSEXP, SEXP all_rhoSEXP, SEXP LR0_allRhoRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type K1R(K1RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K2R(K2RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P0R(P0RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type AR(ARSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type U1R(U1RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wR(wRSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LamdasR(LamdasRSEXP);
    Rcpp::traits::input_parameter< int >::type nminuspx(nminuspxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type all_rho(all_rhoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type LR0_allRhoR(LR0_allRhoRSEXP);
    rcpp_result_gen = Rcpp::wrap(doubleloop_LRT(K1R, K2R, P0R, AR, U1R, wR, LamdasR, nminuspx, all_rho, LR0_allRhoR));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CKLRT_Eigen_C", (DL_FUNC) &_CKLRT_Eigen_C, 1},
    {"_CKLRT_Eigen_C_value", (DL_FUNC) &_CKLRT_Eigen_C_value, 1},
    {"_CKLRT_MatMult_C", (DL_FUNC) &_CKLRT_MatMult_C, 2},
    {"_CKLRT_Sum_C", (DL_FUNC) &_CKLRT_Sum_C, 1},
    {"_CKLRT_ColSum_C", (DL_FUNC) &_CKLRT_ColSum_C, 1},
    {"_CKLRT_MatrixRowMax_C", (DL_FUNC) &_CKLRT_MatrixRowMax_C, 1},
    {"_CKLRT_Elementwisesquare_C", (DL_FUNC) &_CKLRT_Elementwisesquare_C, 1},
    {"_CKLRT_VecMultMat_C", (DL_FUNC) &_CKLRT_VecMultMat_C, 2},
    {"_CKLRT_Vecplus_C", (DL_FUNC) &_CKLRT_Vecplus_C, 2},
    {"_CKLRT_ColSumtwomatrix_C", (DL_FUNC) &_CKLRT_ColSumtwomatrix_C, 2},
    {"_CKLRT_ifelsetest_C", (DL_FUNC) &_CKLRT_ifelsetest_C, 1},
    {"_CKLRT_MatrixPlus_C", (DL_FUNC) &_CKLRT_MatrixPlus_C, 2},
    {"_CKLRT_NumxMatrix_C", (DL_FUNC) &_CKLRT_NumxMatrix_C, 2},
    {"_CKLRT_LR0_fixRho_C", (DL_FUNC) &_CKLRT_LR0_fixRho_C, 5},
    {"_CKLRT_doubleloop", (DL_FUNC) &_CKLRT_doubleloop, 10},
    {"_CKLRT_LR0_fixRho_LRT_C", (DL_FUNC) &_CKLRT_LR0_fixRho_LRT_C, 6},
    {"_CKLRT_doubleloop_LRT", (DL_FUNC) &_CKLRT_doubleloop_LRT, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_CKLRT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
