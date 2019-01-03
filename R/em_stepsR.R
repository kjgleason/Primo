#' E Step
#'
#' DEPRECATED (replaced by faster RcppArmadillo version).
#' Calculate posterior expectations given maximized estimate for \eqn{\pi} (the
#' proportion of SNPs coming from each configuration).
#'
#' @param old_pi vector of configuration proportions, fit through maximization
#' @param n_obs numerical value; number of observations (e.g. SNPs/loci)
#' @param n_pattern numerical value; number of possible patterns
#' (usualy 2^J for J data types)
#' @param q matrix of configurations
#' @param density_0 estimate for the null density function
#' @param density_1 estimate for the alternative density function
#' @param cl cluster object (optional)
#'
#' @return Returns a vector estimating the posterior expectations
#' (i.e. estimated probability of each configuration for each SNP).
#'
#' @export
#'
e_stepR<-function(old_pi, n_obs, n_pattern, q, density_0, density_1, cl=NULL){
  Bmatrix<-matrix(NA,nrow=n_obs,ncol=n_pattern)

  # parallel version
  if(!is.null(cl)){
    Bmatrix <- parallel::parSapply(cl, 1:n_pattern,function(j,q,density_0,density_1,old_pi){
      # calculate log density for current pattern/configuration
      return(as.numeric(log(old_pi[j]) + (log(density_0))%*%(1-q[j,]) +  (log(density_1))%*%q[j,]))
    },q=q, density_0=density_0, density_1=density_1, old_pi=old_pi)
  #sequential version
  }else{
    t_Q <- t(Q)
    Bmatrix <- log(density_0) %*% (1-t_Q) + log(density_1) %*% t_Q
    Bmatrix <- sweep(Bmatrix,2,log(old_pi),"+")
  }

  # subtract minimum, then e^(B)/rowSums(B)
  Bmatrix<-Bmatrix-matrixStats::rowMins(Bmatrix)

  # may need to account for infinite resulting from precision problems
  if(any(Bmatrix>709)){
    Bmatrix <- Rmpfr::mpfr(Bmatrix,precBits=106)
    class(Bmatrix) <- 'mpfrArray'
    Bmatrix<-exp(Bmatrix)
    # can use rowSums directly since columns are recycled
    Bmatrix<-Bmatrix/Rmpfr::rowSums(Bmatrix)
    # convert back to base matrix
    matRows <- nrow(Bmatrix)
    matCols <- ncol(Bmatrix)
    Bmatrix <- matrix(as.numeric(Bmatrix),nrow=matRows,ncol=matCols)
  } else{
    Bmatrix<-exp(Bmatrix)
    # can use rowSums directly since columns are recycled
    Bmatrix<-Bmatrix/matrixStats::rowSums2(Bmatrix)
  }

  return(Bmatrix)
}



#' M Step
#'
#' DEPRECATED (replaced by faster RcppArmadillo version).
#' Estimates the \eqn{\pi} vector that maximizes the posterior expectation function from the E-step.
#'
#' @param old_b matrix of posterior expectations/probabilities ; usually
#' from an expectation (E-)step
#' @param n_obs numerical value; number of observations (e.g. SNPs/loci)
#' @param n_pattern numerical value; number of possible patterns
#' (usually \eqn{2^J} for J data types)
#'
#' @return Returns a vector estimating the proportion of SNPs coming from each configuration.
#'
#' @export
#'
m_stepR<-function(old_b, n_obs, n_pattern){
  # ensure no zero probabilities when calculating pi: newpi<-colSums(old_b)/(n_obs)
  # newpi<-(colSums(old_b)+1)/(n_obs+n_pattern)
  newpi<-colSums(old_b)/n_obs
  return(newpi)
}
