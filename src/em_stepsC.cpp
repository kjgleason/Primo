#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' E Step
//'
//' Calculate posterior expectations given maximized estimate for \eqn{\pi} (the
//' proportion of SNPs coming from each configuration).
//'
//' @param old_pi vector of configuration proportions, fit through maximization
//' (usualy 2^J for J data types)
//' @param Q matrix of configurations
//' @param D_0 estimate for the null density function
//' @param D_1 estimate for the alternative density function
//'
//' @return Returns a vector estimating the posterior expectations
//' (i.e. estimated probability of each configuration for each SNP).
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat e_stepC(const arma::colvec& old_pi, const arma::mat& Q, const arma::mat& D_0, const arma::mat& D_1) {
  // transpose Q to form compatible dimensions for matrix multiplications
  arma::mat t_Q = Q.t();

  // calculate log density under each configuration
  arma::mat Bmat = log(D_0) * (1-t_Q) + log(D_1) * t_Q;
  // factor in proportion of observations estimated to belong to each configuration
  Bmat.each_row() += log(old_pi);

  // substract minimum from each row (i.e. subtract colvec of rowMins from each column)
  Bmat.each_col() -= min(Bmat,1);

  // exponentiate to convert (relative) log density to (relative) density
  Bmat = exp(Bmat);

  // convert to proportion (i.e. divide every column by colvec of rowSums)
  Bmat.each_col() /= sum(Bmat,1);

  return Bmat;
}

//' M Step
//'
//' Estimates the \eqn{\pi} vector that maximizes the posterior expectation function from the E-step.
//'
//' @param old_B matrix of posterior expectations/probabilities ; usually
//' from an expectation (E-)step
//'
//' @return Returns a vector estimating the proportion of SNPs coming from each configuration.
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat m_stepC(const arma::mat& old_B){
  int n_obs = old_B.n_rows, n_pattern = old_B.n_cols;

  arma::mat newpi = (sum(old_B,0)+1)/(n_obs+n_pattern);

  return newpi;
}
