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
  int n_obs = D_0.n_rows, n_pattern = Q.n_rows;

  arma::mat Bmat(n_obs,n_pattern);

  // transpose Q to form compatible dimensions for matrix multiplications
  arma::mat t_Q = Q.t();

  // Consider updating loops to iterator format
  for(int j=0; j < n_pattern; j++){
    Bmat.col(j) = log(old_pi[j]) + log(D_0) * (1-t_Q.col(j)) + log(D_1) * t_Q.col(j);
  }

  for(int i=0; i < n_obs; i++){
    Bmat.row(i)= Bmat.row(i) - Bmat.row(i).min();
  }

  Bmat = exp(Bmat);

  for(int i=0; i < n_obs; i++){
    Bmat.row(i)= Bmat.row(i) / sum(Bmat.row(i));
  }

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
