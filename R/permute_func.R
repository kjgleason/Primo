#' Estimate Posterior Probabilities after Permuting One Column
#'
#' Permute one column and estimate the posterior probability for each configuration
#' for each SNP under the permuted dataset.
#' Utilizes parallel computing, when available.
#'
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param dfs vector of degrees of freedom.
#' @param alt_proportions vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param par_size numerical value; specifies the number of CPUs/cores/processors for parallel computing
#' (0 for sequential processing).
#' @param tol numerical value; specifies tolerance threshold for convergence.
#' @param perm_col numerical value; column number to permute
#' @param density_list (optional) list of densities estimated by \code{estimate_densities()}.
#' @param perm_densities logical value; when true, permutes the previously calculated densities
#' instead of reestimating densities from permuted betas and sds
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
#' \code{alt_proportions} \tab vector of the proportions of test-statistics used in
#' estimating alternative densities\cr
#' \code{tol} \tab numerical value; the tolerance threshold used in determining convergence
#' }
#'
#' @details See documentation for \code{\link{estimate_config}} for additional details
#' regarding the input arguments.
#'
#' @export
#'
permute_once <- function(betas, sds, mafs, dfs, alt_proportions, perm_col, tol=1e-3, par_size=0, density_list=NULL, perm_densities=F){
  # permute the order of rows
  od <- sample(1:nrow(betas))
  # shuffle column using permuted order
  betas[,perm_col] <- betas[od,perm_col]
  sds[,perm_col] <- sds[od,perm_col]
  # match densities of perm_col to permuted order, if densities are provided and perm_densities=T
  if(perm_densities){
    if(is.null(density_list)) stop("Density list cannot be null if perm_densities=T")
    density_list$Tstat_m[,perm_col] <- density_list$Tstat_m[od,perm_col]
    density_list$D0[,perm_col] <- density_list$D0[od,perm_col]
    density_list$D1[,perm_col] <- density_list$D1[od,perm_col]
  } else{
    # calculate new density for permuted data only, if densities are provided
    if(!is.null(density_list)){
      perm_dens <- estimate_densities(betas[,perm_col],sds[,perm_col],mafs,dfs[perm_col],alt_proportions[perm_col])
      density_list$Tstat_m[,perm_col] <- perm_dens$Tstat_m
      density_list$D0[,perm_col] <- perm_dens$D0
      density_list$D1[,perm_col] <- perm_dens$D1
    }
  }
  # estimate configuration for permuted data
  xx <- NULL
  try({xx <- estimate_config(betas, sds, mafs, dfs, alt_proportions=alt_proportions, tol=tol, par_size=par_size, density_list=density_list)}, silent=TRUE)
  return(xx)
}


#' Estimate Posterior Probabilities for Permuted Datasets
#'
#' Permutes several datasets (one column permuted in each)
#' and estimates the posterior probability for each configuration
#' for each SNP under each permuted dataset.
#' Utilizes parallel computing, when available.
#'
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param dfs vector of degrees of freedom.
#' @param alt_proportions vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param tol numerical value; specifies tolerance threshold for convergence.
#' @param par_size numerical value; specifies the number of CPUs/cores/processors for parallel computing
#' (0 for sequential processing).
#' @param perm_par_size numerical value; specifies the number of CPUs/cores/processors for
#' parallel computing of permutations(0 for sequential processing).
#' @param density_list (optional) list of densities estimated by \code{estimate_densities()}.
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
permute_integ <- function(betas, sds, mafs, dfs, alt_proportions, tol=1e-3, par_size=0, perm_par_size=0, density_list=NULL){
  # permute each column
  # parallel version
  if(perm_par_size>0){
    clust <- parallel::makeCluster(perm_par_size)
    # run estimate_config for one permuted dataset on each cluster node
    res <- parallel::parLapply(cl=clust, 1:ncol(betas),
                               function(j,betas,sds,mafs,dfs,alt_proportions,tol,par_size) {
                                 return(permute_once(betas=betas,sds=sds,mafs=mafs,
                                                     dfs=dfs, alt_proportions=alt_proportions,perm_col=j,
                                                     tol=tol,par_size=par_size,density_list=density_list))
                               },betas=betas,sds=sds,mafs=mafs,dfs=dfs,
                               alt_proportions=alt_proportions,tol=tol,par_size=par_size,density_list=density_list)
    parallel::stopCluster(clust)
  # sequential version
  # run estimate_config for each permuted dataset
  } else res <- lapply(1:ncol(betas), function(j)
                        permute_once(betas,sds,mafs,dfs,alt_proportions,perm_col=j,tol,par_size,density_list))

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
#' @param perm_par_size numerical value; specifies the number of CPUs/cores/processors for
#' parallel computing of permutations(0 for sequential processing).
#' @param density_list (optional) list of densities estimated by \code{estimate_densities()}.
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
permute_setup <- function(betas, sds, mafs, dfs, true_res, par_size=0, perm_par_size=0, density_list=NULL){
  # obtain parameters used with true data
  alt_proportions <- true_res$alt_proportions
  tol <- true_res$tol
  if(is.null(density_list)){
    density_list <- list(Tstat_m=true_res$Tstat_m,D0=true_res$D0,D1=true_res$D0)
  }
  # run permutations using same parameters as true data
  res <- permute_integ(betas, sds, mafs, dfs, alt_proportions, tol, par_size, perm_par_size, density_list)
  return(res)
}


#' Estimate Posterior Probabilities for Multiple Permuted Datasets
#'
#' Run multiple permutations of multi-omics datasets.
#' Each set of permutations permutes one column at a time
#' and estimates the posterior probability for each configuration
#' for each SNP under each permuted dataset.
#'
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param dfs vector of degrees of freedom.
#' @param true_res a list; results returned by \code{\link{estimate_config}} when
#' run with true data
#' @param nperm numerical value; number of permutations to run
#' @param par_size numerical value; specifies the number of CPUs/cores/processors for
#' parallel computing (0 for sequential processing).
#' @param perm_par_size numerical value; specifies the number of CPUs/cores/processors for
#' parallel computing of permutations(0 for sequential processing).
#' @param density_list (optional) list of densities estimated by \code{estimate_densities()}.
#'
#' @return A list of lists, where each list holds the matrices of posterior probabilities
#' from \code{\link{estimate_config}} run on a permuted dataset (i.e. a single column).
#'
#' @details See documentation for \code{\link{estimate_config}} for additional details
#' regarding the input arguments.
#'
#' @export
#'
permute_multi <- function(betas, sds, mafs, dfs, true_res, nperm=10, par_size=0, perm_par_size=0, density_list=NULL){
  # obtain parameters used with true data
  alt_proportions <- true_res$alt_proportions
  tol <- true_res$tol
  if(is.null(density_list)){
    density_list <- list(Tstat_m=true_res$Tstat_m,D0=true_res$D0,D1=true_res$D0)
  }

  return_list <- list()
  for(p in 1:nperm){
    # run permutations using same parameters as true data
    res <- permute_integ(betas, sds, mafs, dfs, alt_proportions, tol, par_size, perm_par_size, density_list)
    res <- lapply(res, function(x) x$post_prob)
    return_list[[p]] <- res
  }

  return(return_list)
}
