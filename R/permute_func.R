#' Estimate Posterior Probabilities after Permuting One Column
#'
#' Permute one column and estimate the posterior probability for each configuration
#' for each SNP under the permuted dataset.
#' Utilizes parallel computing, when available.
#'
#' @inheritParams estimate_config
#' @param perm_col numerical value; column number to permute
#'
#' @return A list with the following elements, based on permuted data:
#' \tabular{ll}{
#' \code{post_prob} \tab matrix of posterior probabilities
#' (rows are SNPs; columns are configurations)\cr
#' \code{config_prop} \tab vector of estimated proportion of SNPs
#' belonging to each configuration\cr
#' \code{Tstat_m} \tab matrix of moderated t-statistics\cr
#' \code{D0} \tab matrix of densities calculated under the null distribution\cr
#' \code{D1} \tab matrix of densities calculated under the alternative distribution\cr
#' \code{priors} \tab vector of prior probabilities used in the estimation\cr
#' \code{tol} \tab numerical value; the tolerance threshold used in determining convergence
#' }
#'
#' @details See documentation for \code{\link{estimate_config}} for additional details
#' regarding the input arguments.
#'
#' @export
#'
permute_once <- function(betas, sds, mafs, dfs, priors, perm_col, tol=1e-3, par_size=0){
  # permute the order of rows
  od <- sample(1:nrow(betas))
  # shuffle column using permuted order
  betas[,perm_col] <- betas[od,perm_col]
  sds[,perm_col] <- sds[od,perm_col]
  # estimate configuration for permuted data
  xx <- NULL
  try({xx <- estimate_config(betas, sds, mafs, dfs, priors=priors, tol=tol, par_size=par_size)}, silent=TRUE)
  return(xx)
}


#' Estimate Posterior Probabilities for Permuted Datasets
#'
#' Permutes several datasets (one column permuted in each)
#' and estimates the posterior probability for each configuration
#' for each SNP under each permuted dataset.
#' Utilizes parallel computing, when available.
#'
#' @inheritParams estimate_config
#' @param perm_par_size numerical value; specifies the number of CPUs/cores/processors for
#' parallel computing of permutations(0 for sequential processing).
#'
#' @return A list of lists, where each list holds results from a single run of
#' \code{\link{estimate_config}}. Each run represents results using one permuted dataset.
#' See \code{\link{estimate_config}} for details of each list.
#'
#' @details \code{perm_par_size} differs from \code{par_size} in that permuted datasets
#' are processed in parallel, rather than configurations. Developer recommendation is
#' to parallelize permutations, running configurations sequentially.
#' See documentation for \code{\link{estimate_config}} for additional details
#' regarding the input arguments.
#'
#' @export
#'
permute_integ <- function(betas, sds, mafs, dfs, priors, tol=1e-3, par_size=0, perm_par_size=0){
  # permute each column
  # parallel version
  if(perm_par_size>0){
    clust <- parallel::makeCluster(perm_par_size)
    # run estimate_config for one permuted dataset on each cluster node
    res <- parallel::parLapply(cl=clust, 1:ncol(betas),
                               function(j,betas,sds,mafs,dfs,priors,tol,par_size) {
                                 return(permute_once(betas=betas,sds=sds,mafs=mafs,
                                                     dfs=dfs, priors=priors,perm_col=j,
                                                     tol=tol,par_size=par_size))
                               },betas=betas,sds=sds,mafs=mafs,dfs=dfs,
                               priors=priors,tol=tol,par_size=par_size)
    parallel::stopCluster(clust)
  # sequential version
  # run estimate_config for each permuted dataset
  } else res <- lapply(1:ncol(betas), function(j)
                        permute_once(betas,sds,maf,dfs,priors,perm_col=j,tol,par_size))

  return(res)
}

#' Estimate Posterior Probabilities for Permuted Datasets
#'
#' Wrapper for \code{\link{permute_integ}} that ensures the
#' same parameters are used for permuted datasets as was used
#' when calculating posterior probabilities for the true data.
#' Permutes several datasets (one column permuted in each)
#' and estimates the posterior probability for each configuration
#' for each SNP under each permuted dataset.
#' Utilizes parallel computing, when available.
#'
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param dfs vector of degrees of freedom.
#' @param true_res a list; results returned by \code{\link{estimate_config}} when
#' run with true data
#' @param par_size numerical value; specifies the number of CPUs/cores/processors for
#' parallel computing (0 for sequential processing).
#'
#' @return A list of lists, where each list holds results from a single run
#' of \code{\link{estimate_config}}.Each run represents results using one permuted dataset.
#' See \code{\link{estimate_config}} for details of each list.
#'
#' @details See documentation for \code{\link{estimate_config}} for additional details
#' regarding the input arguments.
#'
#' @export
#'
permute_setup <- function(betas, sds, mafs, dfs, true_res, par_size=0, perm_par_size=0){
  # obtain parameters used with true data
  priors <- true_res$priors
  tol <- true_res$tol
  # run permutations using same parameters as true data
  res <- permute_integ(betas, sds, mafs, dfs, priors, tol, par_size, perm_par_size)
  return(res)
}
