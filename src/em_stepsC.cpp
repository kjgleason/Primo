#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Density matrix calculation
//'
//' Calculate conditional joint densities for each configuration, given marginal
//' null and alternative densities, under the assumption of independence.
//'
//' @param Q matrix of configurations
//' @param D0 matrix of density estimates under the null distribution
//' @param D1 matrix of density estimates under the alternative distributions
//'
//' @return Returns a matrix estimating the conditional joint densities
//' under each configuration
//' (i.e. estimated probability of each configuration for each SNP).
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat calc_Dmatrix(const arma::mat& Q, const arma::mat& D0, const arma::mat& D1) {
  // transpose Q to form compatible dimensions for matrix multiplications
  arma::mat t_Q = Q.t();

  // calculate log density under each configuration
  arma::mat Dmat = log(D0) * (1-t_Q) + log(D1) * t_Q;

  // exponentiate to convert log density to density
  Dmat = exp(Dmat);

  return Dmat;
}

//' E Step
//'
//' Calculate posterior expectations given maximized estimate for \eqn{\pi} (the
//' proportion of SNPs coming from each configuration).
//'
//' @param old_pi vector of configuration proportions, fit through maximization
//' (usualy 2^J for J data types)
//' @param Dmat matrix of conditional joint densities under each configuration
//'
//' @return Returns a matrix estimating the posterior expectations
//' (i.e. estimated probability of each configuration for each SNP).
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat e_step(const arma::rowvec& old_pi, const arma::mat& Dmat) {

  // calculate log density under each configuration
  arma::mat Bmat = log(Dmat);

  // factor in proportion of observations estimated to belong to each configuration
  Bmat.each_row() += log(old_pi);

  // substract maximum from each row (i.e. subtract colvec of rowMaxs from each column)
  Bmat.each_col() -= max(Bmat,1);

  // exponentiate to convert (relative) log density to (relative) density
  Bmat = exp(Bmat);

  // convert to proportion (i.e. divide every column by colvec of rowSums)
  Bmat.each_col() /= sum(Bmat,1);

  return Bmat;
}

//' E Step, using precalculated joint densities
//'
//' Calculate posterior expectations given maximized estimate for \eqn{\pi} (the
//' proportion of SNPs coming from each configuration).
//'
//' @param old_pi vector of configuration proportions, fit through maximization
//' (usualy 2^J for J data types)
//' @param Dmat matrix of conditional joint densities under each configuration
//'
//' @return Returns a matrix of column sums of posterior expectations
//' (allows processing to be performed in chunks).
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat e_step_withColSums(const arma::rowvec& old_pi, const arma::mat& Dmat) {

  // calculate log density under each configuration
  arma::mat Bmat = log(Dmat);
  // factor in proportion of observations estimated to belong to each configuration
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

//' M Step
//'
//' Estimates the \eqn{\pi} vector that maximizes the posterior expectation function from the E-step.
//'
//' @param old_B matrix of posterior expectations/probabilities ; usually
//' from an expectation (E-)step
//'
//' @return Returns a row vector estimating the proportion of SNPs coming from each configuration.
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
//' Calculate posterior expectations given current estimate for \eqn{\pi} (the
//' proportion of SNPs coming from each configuration; E-Step).
//' Re-estimates the \eqn{\pi} vector that maximizes the posterior expectation function from the E-step (i.e. M-step).
//'
//' @param old_pi vector of configuration proportions, fit through maximization
//' (usualy 2^J for J data types)
//' @param Dmat matrix of conditional joint densities under each configuration
//'
//' @return Returns a row vector estimating the proportion of SNPs coming from each configuration.
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat em_iter(const arma::rowvec& old_pi, const arma::mat& Dmat) {
  // E-step

  // calculate log density under each configuration
  arma::mat Bmat = log(Dmat);
  // factor in proportion of observations estimated to belong to each configuration
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
