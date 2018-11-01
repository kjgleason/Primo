#' Estimate Posterior Probabilities after Permuting Densities
#'
#' Permute densities for one or more columns and estimate the posterior probability
#' for each configuration for each SNP under the permuted dataset.
#'
#' @param alt_props vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param perm_col (vector of) numerical value(s); column number(s) to permute
#' @param tol numerical value; specifies tolerance threshold for convergence.
#' @param density_list list of densities estimated by \code{estimate_densities()}.
#' @param config_prop vector of estimated proportion of SNPs
#' belonging to each configuration
#'
#' @return A list with the following elements, based on permuted data:
#' \tabular{ll}{
#' \code{post_prob} \tab matrix of posterior probabilities
#' (rows are SNPs; columns are configurations)\cr
#' \code{config_prop} \tab vector of estimated proportion of SNPs
#' belonging to each configuration\cr
#' \code{Tstat_mod} \tab matrix of moderated t-statistics\cr
#' \code{D0} \tab matrix of densities calculated under the null distribution\cr
#' \code{D1} \tab matrix of densities calculated under the alternative distribution\cr
#' \code{alt_props} \tab vector of the proportions of test-statistics used in
#' estimating alternative densities\cr
#' \code{tol} \tab numerical value; the tolerance threshold used in determining convergence
#' }
#'
#' @details See documentation for \code{\link{estimate_config}} for additional details
#' regarding the input arguments.
#'
#' @export
#'
permute_once_dens <- function(alt_props, perm_col, tol=1e-3, density_list, config_prop){
  m <- nrow(density_list$D0)
  d <- ncol(density_list$D0)

  ## permute densities within requested columns
  for(j in perm_col){
    od <- sample(1:m)
    density_list$Tstat_mod[,j] <- density_list$Tstat_mod[od,j]
    density_list$D0[,j] <- density_list$D0[od,j]
    density_list$D1[,j] <- density_list$D1[od,j]
  }

  # configuration matrix and number of patterns
  Q <- make_qmat(1:d)
  n_pattern <- nrow(Q)

  # Re-calculate posterior probabilities under permutation
  # maximum vector size in R is 2^31-1 (need different calculation if # elements in matrix exceeds)
  if(m*as.double(n_pattern) <= 2^31-1){
    curb<-e_step(config_prop, Q=Q, D0=density_list$D0, D1=density_list$D1)
    # temporary solution to precision problems caused by |t-statistics| > 30
    curb[which(is.na(curb))] <- 0
    return(list(post_prob=curb,config_prop=config_prop,Tstat_mod=density_list$Tstat_mod,
                D0=density_list$D0,D1=density_list$D1,alt_props=alt_props,tol=tol))
  } else{
    curb <- NULL
    exit_e_step <-function(curpi,Q,D0,D1){
      t_Q <- t(Q)
      # consider using `bigmemory` package for following matrix
      curb <- log(D0) %*% (1-t_Q) + log(D1) %*% t_Q
      curb <- sweep(curb,2,log(curpi),"+")
      # subtract minimum, then e^(B)/rowSums(B)
      curb<-curb-matrixStats::rowMins(curb)
      curb<-exp(curb)
      # can use rowSums directly since columns are recycled
      curb<-curb/matrixStats::rowSums2(curb)

      # temporary solution to precision problems caused by |t-statistics| > 30
      ##process rows in chunks (clumsy, but faster than split command)
      n_chunks <- ceiling((m*as.double(n_pattern)) / (2^31-1))
      rowChunks <- list()
      start <- 1
      end <- increment
      for(ch in 1:(n_chunks-1)){
        rowChunks[[ch]] <- start:end
        start <- end+1
        end <- start + (increment-1)
      }
      rowChunks[[n_chunks]] <- start:m

      for(myRow in rowChunks){
        curb[myRow,][which(is.na(curb[myRow,]))] <- 0
      }

      return(curb)
    }
    try(curb<-exit_e_step(curpi=config_prop,Q=Q,D0=density_list$D0,D1=density_list$D1))
    return(list(post_prob=curb,config_prop=config_prop,Tstat_mod=density_list$Tstat_mod,
                D0=density_list$D0,D1=density_list$D1,alt_props=alt_props,tol=tol))
  }
}


#' Estimate Posterior Probabilities after Permuting Column(s) of Statistics
#'
#' Permute one column (beta and sd) and estimate the posterior probability for each configuration
#' for each SNP under the permuted dataset.
#'
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param dfs vector of degrees of freedom.
#' @param alt_props vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param par_size numerical value; specifies the number of CPUs/cores/processors for parallel computing
#' (0 for sequential processing).
#' @param tol numerical value; specifies tolerance threshold for convergence.
#' @param perm_col numerical value; column number to permute
#' @param density_list (optional) list of densities estimated by \code{estimate_densities()}.
#'
#' @return A list with the following elements, based on permuted data:
#' \tabular{ll}{
#' \code{post_prob} \tab matrix of posterior probabilities
#' (rows are SNPs; columns are configurations)\cr
#' \code{config_prop} \tab vector of estimated proportion of SNPs
#' belonging to each configuration\cr
#' \code{Tstat_mod} \tab matrix of moderated t-statistics\cr
#' \code{D0} \tab matrix of densities calculated under the null distribution\cr
#' \code{D1} \tab matrix of densities calculated under the alternative distribution\cr
#' \code{alt_props} \tab vector of the proportions of test-statistics used in
#' estimating alternative densities\cr
#' \code{tol} \tab numerical value; the tolerance threshold used in determining convergence
#' }
#'
#' @details See documentation for \code{\link{estimate_config}} for additional details
#' regarding the input arguments.
#'
#' @export
#'
permute_once_stats <- function(betas, sds, mafs, dfs, alt_props, perm_col, tol=1e-3, par_size=0){
  for(j in perm_col){
    # permute the order of rows
    od <- sample(1:nrow(betas))
    # shuffle column using permuted order
    betas[,j] <- betas[od,j]
    sds[,j] <- sds[od,j]
  }

  # calculate new density for permuted data only, if densities are provided
  if(!is.null(density_list)){
    for(j in perm_col){
      perm_dens <- estimate_densities(betas[,j],sds[,j],mafs,dfs[j],alt_props[j])
      density_list$Tstat_mod[,j] <- perm_dens$Tstat_mod
      density_list$D0[,j] <- perm_dens$D0
      density_list$D1[,j] <- perm_dens$D1
    }
  }

  # estimate configuration for permuted data
  xx <- NULL
  try({xx <- estimate_config(betas, sds, mafs, dfs, alt_props=alt_props, tol=tol, par_size=par_size, density_list=density_list)}, silent=TRUE)
  return(xx)
}


#' Estimate Posterior Probabilities after Permutation
#'
#' Permute column(s) and estimate the posterior probability for each configuration
#' for each SNP under the permuted dataset.
#'
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param dfs vector of degrees of freedom.
#' @param alt_props vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param par_size numerical value; specifies the number of CPUs/cores/processors for parallel computing
#' (0 for sequential processing).
#' @param tol numerical value; specifies tolerance threshold for convergence.
#' @param perm_col numerical value; column number to permute
#' @param density_list list of densities estimated by \code{estimate_densities()}.
#' @param config_prop vector of estimated proportion of SNPs
#' belonging to each configuration
#' @param perm_densities logical value; when true, permutes the previously calculated densities
#' instead of reestimating densities from permuted betas and sds
#'
#' @return A list with the following elements, based on permuted data:
#' \tabular{ll}{
#' \code{post_prob} \tab matrix of posterior probabilities
#' (rows are SNPs; columns are configurations)\cr
#' \code{config_prop} \tab vector of estimated proportion of SNPs
#' belonging to each configuration\cr
#' \code{Tstat_mod} \tab matrix of moderated t-statistics\cr
#' \code{D0} \tab matrix of densities calculated under the null distribution\cr
#' \code{D1} \tab matrix of densities calculated under the alternative distribution\cr
#' \code{alt_props} \tab vector of the proportions of test-statistics used in
#' estimating alternative densities\cr
#' \code{tol} \tab numerical value; the tolerance threshold used in determining convergence
#' }
#'
#' @details See documentation for \code{\link{estimate_config}} for additional details
#' regarding the input arguments.
#'
#' This funciton is a wrapper to call \code{permute_once_dens()} (when \code{perm_densities=TRUE})
#' or \code{permute_once_stats()} (when \code{perm_densities=FALSE}).
#'
#' @export
#'
permute_once <- function(betas=NULL, sds=NULL, mafs=NULL, dfs=NULL, alt_props, perm_col, tol=1e-3, par_size=0,
                         density_list=NULL, config_prop=NULL, perm_densities=T){
  if(perm_densities){
    if(is.null(density_list)) stop("Density list cannot be null if perm_densities=T")
    return(permute_once_dens(alt_props=alt_props, perm_col=perm_col, tol=1e-3, density_list=density_list,config_prop=config_prop))
  } else return(permute_once_stats(betas=betas, sds=sds, mafs=mafs, dfs=dfs, alt_props=alt_props,
                                   perm_col=perm_col, tol=tol, par_size=par_size))
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
#' @param alt_props vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param tol numerical value; specifies tolerance threshold for convergence.
#' @param par_size numerical value; specifies the number of CPUs/cores/processors for parallel computing
#' (0 for sequential processing).
#' @param perm_par_size numerical value; specifies the number of CPUs/cores/processors for
#' parallel computing of permutations(0 for sequential processing).
#' @param density_list (optional) list of densities estimated by \code{estimate_densities()}.
#' @param config_prop vector of estimated proportion of SNPs
#' belonging to each configuration
#' @param perm_densities logical value; when true, permutes the previously calculated densities
#' instead of reestimating densities from permuted betas and sds.
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
permute_integ <- function(betas=NULL, sds=NULL, mafs=NULL, dfs=NULL, alt_props, tol=1e-3, par_size=0, perm_par_size=0,
                          density_list=NULL, config_prop=NULL, perm_densities=T){
  # permute each column
  # parallel version
  if(perm_par_size>0){
    clust <- parallel::makeCluster(perm_par_size)
    parallel::clusterCall(clust, function() {library(primo)})
    # run estimate_config for one permuted dataset on each cluster node
    res <- parallel::parLapply(cl=clust, 1:ncol(betas),
                               function(j,betas,sds,mafs,dfs,alt_props,tol,par_size) {
                                 return(permute_once(betas=betas,sds=sds,mafs=mafs,
                                                     dfs=dfs, alt_props=alt_props,perm_col=j,
                                                     tol=tol,par_size=par_size,density_list=density_list,
                                                     config_prop=config_prop,perm_densities=perm_densities))
                               },betas=betas,sds=sds,mafs=mafs,dfs=dfs,alt_props=alt_props,tol=tol,
                               par_size=par_size,density_list=density_list,config_prop=config_prop,perm_densities=perm_densities)
    parallel::stopCluster(clust)
  # sequential version
  # run estimate_config for each permuted dataset
  } else res <- lapply(1:ncol(betas), function(j)
                        permute_once(betas,sds,mafs,dfs,alt_props,perm_col=j,tol,par_size,density_list,config_prop,perm_densities))

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
#' @param perm_densities logical value; when true, permutes the previously calculated densities
#' instead of reestimating densities from permuted betas and sds.
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
permute_setup <- function(betas=NULL, sds=NULL, mafs=NULL, dfs=NULL, true_res, par_size=0, perm_par_size=0, density_list=NULL, perm_densities=T){
  # obtain parameters used with true data
  alt_props <- true_res$alt_props
  tol <- true_res$tol
  config_prop <- true_res$config_prop
  if(is.null(density_list)){
    density_list <- list(Tstat_mod=true_res$Tstat_mod,D0=true_res$D0,D1=true_res$D0)
  }
  # run permutations using same parameters as true data
  res <- permute_integ(betas, sds, mafs, dfs, alt_props, tol, par_size, perm_par_size, density_list, config_prop, perm_densities=T)
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
#' @param perm_densities logical value; when true, permutes the previously calculated densities
#' instead of reestimating densities from permuted betas and sds.
#'
#' @return A list of lists, where each list holds the matrices of posterior probabilities
#' from \code{\link{estimate_config}} run on a permuted dataset (i.e. a single column).
#'
#' @details See documentation for \code{\link{estimate_config}} for additional details
#' regarding the input arguments.
#'
#' @export
#'
permute_multi <- function(betas=NULL, sds=NULL, mafs=NULL, dfs=NULL, true_res, nperm=10, par_size=0, perm_par_size=0, density_list=NULL, perm_densities=T){
  # obtain parameters used with true data
  alt_props <- true_res$alt_props
  tol <- true_res$tol
  config_prop <- true_res$config_prop
  if(is.null(density_list)){
    density_list <- list(Tstat_mod=true_res$Tstat_mod,D0=true_res$D0,D1=true_res$D0)
  }

  return_list <- list()
  for(p in 1:nperm){
    # run permutations using same parameters as true data
    res <- permute_integ(betas, sds, mafs, dfs, alt_props, tol, par_size, perm_par_size, density_list, config_prop, perm_densities)
    res <- lapply(res, function(x) x$post_prob)
    return_list[[p]] <- res
  }

  return(return_list)
}
