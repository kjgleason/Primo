#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Density matrix calculation
//'
//' Calculate conditional joint densities for each association pattern, given marginal
//' null and alternative densities, under the assumption of independence.
//'
//' @param Q matrix of association patterns.
//' @param D0 matrix of density estimates under the null distribution.
//' @param D1 matrix of density estimates under the alternative distributions.
//'
//' @return A matrix estimating the conditional joint densities
//' under each association pattern for each observation (e.g SNP).
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat calc_Dmatrix(const arma::mat& Q, const arma::mat& D0, const arma::mat& D1) {
  // transpose Q to form compatible dimensions for matrix multiplications
  arma::mat t_Q = Q.t();

  // calculate log density under each association pattern
  arma::mat Dmat = log(D0) * (1-t_Q) + log(D1) * t_Q;

  // exponentiate to convert log density to density
  Dmat = exp(Dmat);

  return Dmat;
}

//' E-step
//'
//' Calculate posterior expectations given maximized estimate for \eqn{\pi}, the
//' proportion of observations (e.g. SNPs) coming from each association pattern.
//'
//' @param old_pi vector of association pattern proportions, fit through maximization.
//' @param Dmat matrix of conditional joint densities under each association pattern.
//'
//' @return A matrix estimating the posterior expectations
//' (i.e. estimated proportion of observations belonging to each association pattern).
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat e_step(const arma::rowvec& old_pi, const arma::mat& Dmat) {

  // calculate log density under each association pattern
  arma::mat Bmat = log(Dmat);

  // factor in proportion of observations estimated to belong to each association patterns
  Bmat.each_row() += log(old_pi);

  // substract maximum from each row (i.e. subtract colvec of rowMaxs from each column)
  Bmat.each_col() -= max(Bmat,1);

  // exponentiate to convert (relative) log density to (relative) density
  Bmat = exp(Bmat);

  // convert to proportion (i.e. divide every column by colvec of rowSums)
  Bmat.each_col() /= sum(Bmat,1);

  return Bmat;
}

//' E-step, with column sums of posterior probabilities
//'
//' Calculate posterior expectations given maximized estimate for \eqn{\pi} (the
//' proportion of observations coming from each association pattern). Return column sums of the posterior
//' expectations. This allows the M-step to be completed in chunks (dividing a running sum
//' by the total number of rows at the end to obtain the mean), alleviating potential memory issues.
//'
//' @inheritParams e_step
//'
//' @return A row vector of column sums of posterior expectations
//' (allows M-step to be performed in chunks).
//'
//' @details
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat e_step_withColSums(const arma::rowvec& old_pi, const arma::mat& Dmat) {

  // calculate log density under each association pattern
  arma::mat Bmat = log(Dmat);
  // factor in proportion of observations estimated to belong to each association pattern
  Bmat.each_row() += log(old_pi);

  // substract maximum from each row (i.e. subtract colvec of rowMaxs from each column)
  Bmat.each_col() -= max(Bmat,1);

  // exponentiate to convert (relative) log density to (relative) density
  Bmat = exp(Bmat);

  // convert to proportion (i.e. divide every column by colvec of rowSums)
  Bmat.each_col() /= sum(Bmat,1);

  // return colSums of posterior probabilities
  arma::mat Bmat_colsums = sum(Bmat,0);

  return Bmat_colsums;
}

//' M-step
//'
//' Estimates the \eqn{\pi} vector that maximizes the posterior expectation function from the E-step.
//'
//' @param old_B matrix of posterior expectations/probabilities (e.g. from an E-step).
//'
//' @return A row vector estimating the proportion of observations coming from each association pattern.
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat m_step(const arma::mat& old_B){
  int n_obs = old_B.n_rows;

  arma::mat newpi = sum(old_B,0)/n_obs;

  return newpi;
}

//' EM Iteration
//'
//' Complete one iteration of the EM-algorithm (combines E-step and M-step).
//' Calculate posterior expectations given current estimate for \eqn{\pi},the
//' proportion of observations (e.g. SNPs) coming from each association pattern
//' (i.e. E-Step). Re-estimate the \eqn{\pi} vector that maximizes the posterior
//' expectations (i.e. M-step).
//'
//' @inheritParams e_step
//'
//' @inherit m_step return
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat em_iter(const arma::rowvec& old_pi, const arma::mat& Dmat) {
  // E-step

  // calculate log density under each association pattern
  arma::mat Bmat = log(Dmat);
  // factor in proportion of observations estimated to belong to each association pattern
  Bmat.each_row() += log(old_pi);

  // substract maximum from each row (i.e. subtract colvec of rowMins from each column)
  Bmat.each_col() -= max(Bmat,1);

  // exponentiate to convert (relative) log density to (relative) density
  Bmat = exp(Bmat);

  // convert to proportion (i.e. divide every column by colvec of rowSums)
  Bmat.each_col() /= sum(Bmat,1);


  // M-step
  int n_obs = Bmat.n_rows;
  // Calculate pi vector that maximizes posterior expectation
  arma::mat newpi = sum(Bmat,0)/n_obs;

  return newpi;
}
