#' Estimate posterior probabilities of association patterns, using t-statistics.
#'
#' For each SNP, estimates the posterior probability for each configuration.
#' Utilizes parallel computing, when available.
#'
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param dfs vector or matrix of degrees of freedom.
#' @param alt_props vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param mafs vector or matrix of minor allele frequencies (MAFs).
#' @param Gamma correlation matrix.
#' @param tol numeric value; specifies tolerance threshold for convergence.
#' @param par_size numeric value; specifies the number of CPUs/cores/processors for
#' parallel computing (1 for sequential processing).
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{post_prob} \tab matrix of posterior probabilities
#' (each column corresponds to an association pattern).\cr
#' \code{pis} \tab vector of estimated proportion of observations
#' belonging to each association pattern.\cr
#' \code{D_mat} \tab matrix of densities under each association pattern.\cr
#' \code{Gamma} \tab correlation matrix.\cr
#' \code{Tstat_mod} \tab matrix of moderated t-statistics.\cr
#' \code{V_mat} \tab matrix of scaling factors under the alternative distribution.\cr
#' \code{mdf_sd_mat} \tab matrix of standard deviation adjustment according to
#'  moderated degrees of freedom: df/(df-2).\cr
#' \code{prior_df} \tab vector of the prior degrees of freedom for each marginal distribution.\cr
#' \code{prior_var} \tab vector of the prior variance estimators for each marginaldistribution.\cr
#' \code{unscaled_var} \tab vector of the unscaled variance priors on non-zero coefficients
#' for each marginal distribution.
#' }
#'
#' The main element of interest for inference is the posterior probabilities matrix, \code{post_prob}.
#' The estimated proportion of observations belonging to each association pattern, \code{pis}, may
#' also be of interest. The remaining elements are returned primarily for use by other functions --
#' such as those conducting conditional association analysis.
#'
#' @details The following are additional details describing the input arguments
#'  (for \eqn{m} SNPs/observations measured in \eqn{d} studies):
#' \tabular{ll}{
#' \code{betas} \tab  \eqn{m} x \eqn{d} matrix.\cr
#' \code{sds} \tab \eqn{m} x \eqn{d} matrix.\cr
#' \code{dfs} \tab vector of length \eqn{d} or an \eqn{m} x \eqn{d} matrix.\cr
#' \code{alt_props} \tab vector of length \eqn{d}.\cr
#' \code{mafs} \tab vector of length \eqn{m} or an \eqn{m} x \eqn{d} matrix.\cr
#'  \tab If \code{NULL}, standard errors will not be adjusted for MAF.\cr
#' \code{Gamma} \tab  \eqn{d} x \eqn{d} matrix.\cr
#'  \tab If \code{NULL}, will be estimated using observations where all \eqn{|t| < 5}.\cr
#' }
#'
#' @export
#'
Primo_tstat <- function(betas, sds,  dfs, alt_props, mafs=NULL, Gamma=NULL, tol=0.001,par_size=1){
  m <- nrow(betas)
  d <- ncol(betas)

  if (is.null(Gamma)) {
    tt <- betas/sds
    tt[abs(tt)>=5] <- NA
    Gamma<- cor(tt,use="complete")
  }

  ## estimate marginal density functions in limma framework
  density_list <- lapply(1:d, function(j){
    if(is.matrix(mafs)) mafs <- mafs[,j]
    if(is.matrix(dfs)) {
      dfs <- dfs[,j]
    } else dfs <- rep(dfs[j],m)
    primo::estimate_densities_modT(betas=betas[,j],sds=sds[,j],mafs=mafs,df=dfs,alt_prop=alt_props[j])
  } )

  ## extract parameters for pattern-specific density estimation
  Tstat_mod <- do.call("cbind", lapply(density_list, function(x) x$Tstat_mod))
  mdfs <- do.call("cbind", lapply(density_list, function(x) x$df_mod))
  V <- do.call("cbind", lapply(density_list, function(x) x$scaler))
  ## extract parameters from marginal densities (can be used to estimate variables for new observations)
  prior_df <- sapply(density_list, function(x) x$prior_df)
  prior_var <- sapply(density_list, function(x) x$prior_var)
  unscaled_var <- sapply(density_list, function(x) x$unscaled_var)

  ## create sd matrix from moderated t-statistic dfs
  mdf_sd_mat <- sqrt(mdfs/(mdfs-2))

  ## computation of D_mat (densities under each pattern)
  Q<-make_qmat(1:d)
  D_mat_func <- function(k){
    m=nrow(V)
    d=ncol(V)
    v_k <-   V %*%diag(Q[k,]) + matrix(1, nrow=nrow(V),ncol=ncol(V))%*%diag( 1-Q[k,])
    v_k <- v_k * mdf_sd_mat
    Z <- Tstat_mod/v_k
    dets <- exp(rowSums(log(v_k)))
    return(mvtnorm::dmvnorm(Z, mean=rep(0,d),sigma=Gamma)/dets)
  }

  ## parallel version
  if(par_size > 1){
    cl <- parallel::makeCluster(par_size)
    parallel::clusterEvalQ(cl, library(mvtnorm))
    parallel::clusterExport(cl=cl, varlist=list("D_mat_func","Tstat_mod","V","mdf_sd_mat","Q","Gamma"),envir=environment())
    D_mat_byk <- parallel::parLapply(cl, 1:2^d, function(k) D_mat_func(k=k))
    parallel::stopCluster(cl)
  } else{
    ## sequential version
    D_mat_byk <- lapply(1:2^d, function(k) D_mat_func(k=k))
  }

  ## bind densities together into matrix
  D_mat <- do.call("cbind",D_mat_byk)
  rm(D_mat_byk)

  ## resume algorithm now that D_mat has been calculated either in parallel or sequentially
  n_pattern <- 2^d
  curpi<- c(0.80, rep((1-0.80)/(n_pattern-1),n_pattern-1))
  diff<-1
  numiters<-1
  itermat<-curpi

  ## Rcpp has maximum vector size of 2^31-1
  ## process large matrices in chunks
  if(m*as.double(n_pattern) <= 2^31-1){
    while(diff>tol){
      start_time <- Sys.time()
      cat("\nIteration:",numiters)
      numiters<-numiters+1
      ## one iteration of EM
      curpi <- primo::em_iter(curpi,D_mat)
      itermat<-rbind(itermat,curpi)
      diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
      Sys.time() - start_time
    }

    PP<- primo::e_step(curpi, Dmat=D_mat)

  } else{

    num_chunks <- ceiling((m*as.double(n_pattern))/(2^31-1))
    Drow_chunks <- split(1:m, ceiling(seq_along(1:m)/(m/num_chunks)))

    while(diff>tol){
      start_time <- Sys.time()
      cat("\nIteration:",numiters)
      numiters<-numiters+1
      # e-step, in chunks
      curb_colsums <- 0
      for(ch in 1:length(Drow_chunks)){
        D_rows <-  Drow_chunks[[ch]]
        # EM on current chunk
        curb_colsums_temp<-primo::e_step_withColSums(curpi, Dmat=D_mat[D_rows,])
        curb_colsums <- curb_colsums + curb_colsums_temp
      }
      curpi<-curb_colsums/m
      itermat<-rbind(itermat,curpi)
      diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
      Sys.time() - start_time
    }

    ## obtain posterior probabilities
    PP <- sweep(log(D_mat),2,log(curpi),"+")
    # subtract maximum, then e^(B)/rowSums(B)
    PP<-PP-matrixStats::rowMaxs(PP)
    PP<-exp(PP)
    # use rowSums directly since columns are recycled
    PP<-PP/rowSums(PP)

  }

  return(list(post_prob=PP, pis=curpi, D_mat=D_mat, Gamma=Gamma, Tstat_mod=Tstat_mod, V_mat = V, mdf_sd_mat = mdf_sd_mat,
              prior_df=prior_df,prior_var=prior_var,unscaled_var=unscaled_var))
}


#' Estimate posterior probabilities of association patterns, using P-values.
#'
#' For each SNP, estimates the posterior probability for each configuration.
#' Utilizes parallel computing, when available.
#'
#' @param pvals matrix of \eqn{P}-values from test statistics.
#' @param alt_props vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param Gamma correlation matrix.
#' @param tol numeric value; specifies tolerance threshold for convergence.
#' @param par_size numeric value; specifies the number of CPUs/cores/processors for
#' parallel computing (1 for sequential processing).
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{post_prob} \tab matrix of posterior probabilities
#' (rows are SNPs; columns are association patterns).\cr
#' \code{pis} \tab vector of estimated proportion of SNPs
#' belonging to each association pattern.\cr
#' \code{D_mat} \tab matrix of densities under each association pattern.\cr
#' \code{Gamma} \tab correlation matrix.\cr
#' \code{chi_mix} \tab matrix of \eqn{-2}log(\eqn{P})-values.\cr
#' \code{A} \tab vector of scaling factors under the alternative distributions.\cr
#' \code{df_alt} \tab vector of degrees of freedom approximated for the alternative distributions.\cr
#' }
#'
#' @details The following are additional details describing the input arguments
#'  (for \eqn{m} SNPs/observations measured in \eqn{d} studies):
#' \tabular{ll}{
#' \code{pvals} \tab  \eqn{m} x \eqn{d} matrix.\cr
#' \code{alt_props} \tab vector of length \eqn{d}.\cr
#' \code{Gamma} \tab  \eqn{d} x \eqn{d} matrix.
#'  If \code{NULL}, will be estimated using observations where all \eqn{p < 5.7e-7}.\cr
#' }
#'
#' @export
#'
Primo_pval <- function(pvals, alt_props, Gamma=NULL, tol=0.001, par_size=1){
  m <- nrow(pvals)
  d <- ncol(pvals)

  chi_mix<-(-2)*log(pvals)

  if (is.null(Gamma)) {
    xx <- chi_mix
    xx[pvals < 5.7e-7] <- NA ## 5.7 * 10^-7 is approx p-val for abs(t)=5
    Gamma<- cor(xx,use="complete")
  }

  # A <- NULL
  # df_alt <- NULL
  #
  # ## consider parallelizing this step for large d
  # ## estimate scaling factor and degrees of freedom for each alternative distribution
  # for(j in 1:d){
  #   optim_dat <- list(chi_mix=sort(chi_mix[,j], decreasing=T), alt_props=alt_props[j])
  #   init1 <- 2
  #   init2 <- 3
  #   ## run global optimization to identify correct "neighborhood" of optimum
  #   global_res <- nloptr::nloptr(x0=c(init1,init2),eval_f=chiMix_pDiff,lb=c(1,2),ub=c(100,100),
  #                                opts=list(algorithm="NLOPT_GN_DIRECT",maxeval=500),
  #                                data=optim_dat, sorted=T)
  #   init1 <- global_res$solution[1]
  #   init2 <- global_res$solution[2]
  #   ## refine optimum using a local optimization algorithm
  #   optim_res <- nloptr::nloptr(x0=c(init1,init2),eval_f=chiMix_pDiff,lb=c(1,2),
  #                               opts=list(algorithm="NLOPT_LN_COBYLA",maxeval=500),
  #                               data=optim_dat, sorted=T)
  #
  #   A <- c(A,optim_res$solution[1])           ## store scale parameter
  #   df_alt <- c(df_alt,optim_res$solution[2]) ## store degrees of freedom
  # }

  ## estimate marginal density functions in limma framework
  density_list <- lapply(1:d, function(j){
    primo::estimate_densities_pval(pvals=pvals[,j],alt_prop=alt_props[j])
  } )

  ## extract parameters from marginal densities
  A <- sapply(density_list, function(x) x$A)
  df_alt <- sapply(density_list, function(x) x$df_alt)

  ## computation of D_mat (densities under each pattern)
  Q<-make_qmat(1:d)
  D_mat_func_p <- function(k){
    A_k <- A*Q[k,] + 1*(1-Q[k,])
    df_k <- df_alt*Q[k,] + 2*(1-Q[k,])

    rate_k <- 1/(2*A_k)
    shape_k <- df_k/2

    return(lcmix::dmvgamma(chi_mix, shape=shape_k, rate=rate_k, corr=Gamma))
  }

  ## parallel version
  if(par_size > 1){
    cl <- parallel::makeCluster(par_size)
    parallel::clusterEvalQ(cl, library(lcmix))
    parallel::clusterExport(cl=cl, varlist=list("D_mat_func_p","chi_mix","A","df_alt","Q","Gamma"),envir=environment())
    D_mat_byk <- parallel::parLapply(cl, 1:2^d, function(k) D_mat_func_p(k=k))
    parallel::stopCluster(cl)
  }else{
    ## sequential version
    D_mat_byk <- lapply(1:2^d, function(k) D_mat_func_p(k=k))
  }

  D_mat <- do.call("cbind",D_mat_byk)
  rm(D_mat_byk)

  ## TO DO: develop more sophisticated handling of "zero" densities
  D_mat[which(D_mat==0)] <- min(D_mat[which(D_mat != 0)])

  n_pattern <- 2^d
  curpi<- c(0.80, rep((1-0.80)/(2^d-1),2^d-1))
  diff<-1
  numiters<-1
  itermat<-curpi

  ## Rcpp has maximum vector size of 2^31-1
  ## process large matrices in chunks
  if(m*as.double(n_pattern) <= 2^31-1){
    while(diff>tol){
      start_time <- Sys.time()
      cat("\nIteration:",numiters)
      numiters<-numiters+1
      ## one iteration of EM
      curpi <- primo::em_iter(curpi,D_mat)
      itermat<-rbind(itermat,curpi)
      diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
      Sys.time() - start_time
    }

    PP<- primo::e_step(curpi, Dmat=D_mat)

  } else{

    num_chunks <- ceiling((m*as.double(n_pattern))/(2^31-1))
    Drow_chunks <- split(1:m, ceiling(seq_along(1:m)/(m/num_chunks)))

    while(diff>tol){
      start_time <- Sys.time()
      cat("\nIteration:",numiters)
      numiters<-numiters+1
      # e-step, in chunks
      curb_colsums <- 0
      for(ch in 1:length(Drow_chunks)){
        D_rows <-  Drow_chunks[[ch]]
        # EM on current chunk
        curb_colsums_temp<-primo::e_step_withColSums(curpi, Dmat=D_mat[D_rows,])
        curb_colsums <- curb_colsums + curb_colsums_temp
      }
      curpi<-curb_colsums/m
      itermat<-rbind(itermat,curpi)
      diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
      Sys.time() - start_time
    }

    ## obtain posterior probabilities
    PP <- sweep(log(D_mat),2,log(curpi),"+")
    # subtract maximum, then e^(B)/rowSums(B)
    PP<-PP-matrixStats::rowMaxs(PP)
    PP<-exp(PP)
    # use rowSums directly since columns are recycled
    PP<-PP/rowSums(PP)

  }

  return(list(postprob=PP, pis=curpi, D_mat=D_mat, Gamma=Gamma, chi_mix=chi_mix, A=A, df_alt=df_alt))
}


#' Estimate posterior probabilities of association patterns, using moderated t-statistics.
#'
#' This version of the function uses moderated \eqn{t}-statistics and parameters
#' previously calculated under the limma framework. It is useful for cases where
#' the same statistic from one study (e.g. gene-SNP pair) may be mapped to
#' multiple statistics from another study (e.g. multiple gene-CpG pairs).
#' For each SNP, estimates the posterior probability for each configuration.
#' Utilizes parallel computing, when available.
#'
#' @param Tstat_mod matrix of moderated t-statistics.
#' @param mdfs matrix of moderated degrees of freedom.
#' @param V matrix of scaling factors.
#' @param Gamma correlation matrix.
#' @param tol numerical value; specifies tolerance threshold for convergence.
#' @param par_size numerical value; specifies the number of CPUs/cores/processors for
#' parallel computing (1 for sequential processing).
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{post_prob} \tab matrix of posterior probabilities
#' (rows are SNPs; columns are association patterns).\cr
#' \code{pis} \tab vector of estimated proportion of SNPs
#' belonging to each association pattern.\cr
#' \code{D_mat} \tab matrix of densities under each association pattern.\cr
#' \code{Gamma} \tab correlation matrix.\cr
#' \code{Tstat_mod} \tab matrix of moderated t-statistics.\cr
#' \code{V_mat} \tab matrix of scaling factors under the alternative distribution.\cr
#' \code{mdf_sd_mat} \tab matrix of standard deviation adjustment according to
#'  moderated degrees of freedom: df/(df-2).\cr
#' }
#'
#' @details The following are additional details describing the input arguments
#'  (for \eqn{m} SNPs/observations measured in \eqn{d} studies):
#' \tabular{ll}{
#' \code{Tstat_mod} \tab  \eqn{m} x \eqn{d} matrix.\cr
#' \code{mdfs} \tab \eqn{m} x \eqn{d} matrix.\cr
#' \code{V} \tab \eqn{m} x \eqn{d} matrix.\cr
#' \code{Gamma} \tab  \eqn{d} x \eqn{d} matrix.\cr
#' }
#'
#' @export
#'
Primo_ModT <- function(Tstat_mod, mdfs, V_mat, Gamma, tol=0.001,par_size=1){
  m <- nrow(Tstat_mod)
  d <- ncol(Tstat_mod)

  ## when degrees of freedom are the same for 1 phenotype across observations, rbind sd based on mdf into matrix format
  if(!is.matrix(mdfs)){
    mdf_sd_mat <- matrix(rep(sqrt(mdfs/(mdfs-2)),each=m),ncol=d)
  } else{
    mdf_sd_mat <- sqrt(mdfs/(mdfs-2))
  }

  ## parallel computation of D_mat(densities under each pattern)
  Q<-make_qmat(1:d)
  D_mat_func <- function(k){
    m=nrow(V_mat)
    d=ncol(V_mat)
    v_k <-   V_mat %*%diag(Q[k,]) + matrix(1, nrow=nrow(V_mat),ncol=ncol(V_mat))%*%diag( 1-Q[k,])
    v_k <- v_k * mdf_sd_mat
    Z <- Tstat_mod/v_k
    dets <- exp(rowSums(log(v_k)))*sqrt(det(Gamma))
    return(mvtnorm::dmvnorm(Z, mean=rep(0,d),sigma=Gamma)/dets)
  }

  ## parallel version
  if(par_size > 1){
    cl <- parallel::makeCluster(par_size)
    parallel::clusterEvalQ(cl, library(mvtnorm))
    parallel::clusterExport(cl=cl, varlist=list("D_mat_func","Tstat_mod","V_mat","mdf_sd_mat","Q","Gamma"),envir=environment())
    D_mat_byk <- parallel::parLapply(cl, 1:2^d, function(k) D_mat_func(k=k))
    parallel::stopCluster(cl)
  } else{
    ## sequential version
    D_mat_byk <- lapply(1:2^d, function(k) D_mat_func(k=k))
  }

  D_mat <- do.call("cbind",D_mat_byk)
  rm(D_mat_byk)

  ## resume algorithm now that D_mat has been calculated either in parallel or sequentially
  n_pattern <- 2^d
  curpi<- c(0.80, rep((1-0.80)/(n_pattern-1),n_pattern-1))
  diff<-1
  numiters<-1
  itermat<-curpi

  ## Rcpp has maximum vector size of 2^31-1
  ## for now, process large matrices in chunks (may eventually utilize bigmemory package)
  if(m*as.double(n_pattern) <= 2^31-1){
    while(diff>tol){
      start_time <- Sys.time()
      cat("\nIteration:",numiters)
      numiters<-numiters+1
      curpi <- primo::em_iter(curpi,D_mat)
      itermat<-rbind(itermat,curpi)
      diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
      Sys.time() - start_time
    }

    PP<- primo::e_step(curpi, Dmat=D_mat)

  } else{

    num_chunks <- ceiling((m*as.double(n_pattern))/(2^31-1))
    Drow_chunks <- split(1:m, ceiling(seq_along(1:m)/(m/num_chunks)))

    while(diff>tol){
      start_time <- Sys.time()
      cat("\nIteration:",numiters)
      numiters<-numiters+1
      # e-step, in chunks
      curb_colsums <- 0
      for(ch in 1:length(Drow_chunks)){
        D_rows <-  Drow_chunks[[ch]]
        curb_colsums_temp<-primo::e_step_withColSums(curpi, Dmat=D_mat[D_rows,])
        curb_colsums <- curb_colsums + curb_colsums_temp
      }
      curpi<-curb_colsums/m
      itermat<-rbind(itermat,curpi)
      diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
      Sys.time() - start_time
    }

    ## obtain posterior probabilities
    PP <- sweep(log(D_mat),2,log(curpi),"+")
    # subtract maximum, then e^(B)/rowSums(B)
    PP<-PP-matrixStats::rowMaxs(PP)
    PP<-exp(PP)
    # can use rowSums directly since columns are recycled
    PP<-PP/matrixStats::rowSums2(PP)

  }

  return(list(post_prob=PP, pis=curpi, D_mat=D_mat, Gamma=Gamma, Tstat_mod=Tstat_mod, V_mat = V_mat, mdf_sd_mat = mdf_sd_mat))
}


#' Estimate posterior probabilities of association patterns.
#'
#' For each observation (e.g. SNP), estimate the posterior probability for each association pattern.
#' This function serves as a wrapper to call either the t-statistic or p-value version of Primo.
#' Utilizes parallel computing, when available.
#'
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param dfs vector or matrix of degrees of freedom.
#' @param pvals matrix of p-values.
#' @param alt_props vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param mafs vector or matrix of minor allele frequencies (MAFs).
#' @param Gamma correlation matrix.
#' @param tol numeric value; specifies tolerance threshold for convergence.
#' @param par_size numeric value; specifies the number of CPUs/cores/processors for
#' parallel computing (1 for sequential processing).
#' @param use_method character string denoting which method to use.
#' Must be one of "tstat" or "pval".
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{post_prob} \tab matrix of posterior probabilities
#' (each column corresponds to an association pattern).\cr
#' \code{pis} \tab vector of estimated proportion of observations
#' belonging to each association pattern.\cr
#' \code{D_mat} \tab matrix of densities under each association pattern.\cr
#' \code{Gamma} \tab correlation matrix.\cr
#' }
#'
#' \itemize{
#' \item If \code{use_method}="tstat", the list will additionally contain:
#' \tabular{ll}{
#' \code{Tstat_mod} \tab matrix of moderated t-statistics.\cr
#' \code{V_mat} \tab matrix of scaling factors under the alternative distribution.\cr
#' \code{mdf_sd_mat} \tab matrix of standard deviation adjustment according to
#'  moderated degrees of freedom: df/(df-2).\cr
#' \code{prior_df} \tab vector of the prior degrees of freedom for each marginal distribution.\cr
#' \code{prior_var} \tab vector of the prior variance estimators for each marginaldistribution.\cr
#' \code{unscaled_var} \tab vector of the unscaled variance priors on non-zero coefficients
#' for each marginal distribution.
#' }
#'
#' \item If \code{use_method}="pval", the list will additionally contain:
#' \tabular{ll}{
#' \code{chi_mix} \tab matrix of \eqn{-2}log(\eqn{P})-values.\cr
#' \code{A} \tab vector of scaling factors under the alternative distributions.\cr
#' \code{df_alt} \tab vector of degrees of freedom approximated for the alternative distributions.\cr
#'  }
#'  }
#'
#' The primary element of interest for inference is the posterior probabilities matrix, \code{post_prob}.
#' The estimated proportion of observations belonging to each association pattern, \code{pis}, may
#' also be of interest.
#' The remaining elements are returned primarily for use by other functions --
#' such as those conducting conditional association analysis.
#'
#' @details The following are additional details describing the input arguments
#'  (for \eqn{m} SNPs/observations measured in \eqn{d} studies):
#' \tabular{ll}{
#' \code{betas} \tab  \eqn{m} x \eqn{d} matrix.\cr
#' \code{sds} \tab \eqn{m} x \eqn{d} matrix.\cr
#' \code{dfs} \tab vector of length \eqn{d} or an \eqn{m} x \eqn{d} matrix.\cr
#' \code{pvals} \tab \eqn{m} x \eqn{d} matrix.\cr
#' \code{alt_props} \tab vector of length \eqn{d}.\cr
#' \code{mafs} \tab vector of length \eqn{m} or an \eqn{m} x \eqn{d} matrix.\cr
#'  \tab If \code{NULL}, standard errors will not be adjusted for MAF.\cr
#' \code{Gamma} \tab  \eqn{d} x \eqn{d} matrix.\cr
#'  \tab If \code{NULL}, will be estimated using observations where all \eqn{|t| < 5}
#'  (or all \eqn{P < 5.7e-7}).\cr
#' }
#'
#' @export
#'
Primo <- function(betas=NULL, sds=NULL, dfs=NULL, pvals=NULL, alt_props, mafs=NULL, Gamma=NULL, tol=1e-3, par_size=1, use_method="tstat"){

  if(use_method=="tstat"){

    if(is.null(betas) | is.null(sds) | is.null(dfs)) stop("When use_method='tstat', arguments betas, sds and dfs must all be non-NULL.")

    return(Primo_tstat(betas, sds,  dfs, alt_props, mafs, Gamma, tol,par_size))

  } else if(use_method=="pval"){

    if(is.null(pvals)) stop("When use_method='pval', argument pvals cannot be NULL.")

    return(Primo_pval(pvals, alt_props, Gamma, tol))

  } else{

    stop("use_method must be either 'tstat' or 'pvalue'.")

  }

}
