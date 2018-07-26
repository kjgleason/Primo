// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// e_step
arma::mat e_step(const arma::rowvec& old_pi, const arma::mat& Q, const arma::mat& D_0, const arma::mat& D_1);
RcppExport SEXP _primo_e_step(SEXP old_piSEXP, SEXP QSEXP, SEXP D_0SEXP, SEXP D_1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type old_pi(old_piSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D_0(D_0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D_1(D_1SEXP);
    rcpp_result_gen = Rcpp::wrap(e_step(old_pi, Q, D_0, D_1));
    return rcpp_result_gen;
END_RCPP
}
// m_step
arma::mat m_step(const arma::mat& old_B);
RcppExport SEXP _primo_m_step(SEXP old_BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type old_B(old_BSEXP);
    rcpp_result_gen = Rcpp::wrap(m_step(old_B));
    return rcpp_result_gen;
END_RCPP
}
// findDiffGit
arma::mat findDiffGit(const arma::rowvec& old_pi, const arma::mat& Q, const arma::mat& D_0, const arma::mat& D_1);
RcppExport SEXP _primo_findDiffGit(SEXP old_piSEXP, SEXP QSEXP, SEXP D_0SEXP, SEXP D_1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type old_pi(old_piSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D_0(D_0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D_1(D_1SEXP);
    rcpp_result_gen = Rcpp::wrap(findDiffGit(old_pi, Q, D_0, D_1));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_primo_e_step", (DL_FUNC) &_primo_e_step, 4},
    {"_primo_m_step", (DL_FUNC) &_primo_m_step, 1},
    {"_primo_findDiffGit", (DL_FUNC) &_primo_findDiffGit, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_primo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
