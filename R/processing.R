#' Subset a Primo object.
#'
#' Subset results from Primo output based on a vector of indices.
#'
#' @param Primo_obj list of results returned by Primo (from the function
#' \code{\link{Primo_tstat}}, \code{\link{Primo_pval}}, \code{\link{Primo_modT}}
#' or \code{\link{Primo_chiMix}}).
#' @param idx integer vector of the indices to subset Primo results.
#'
#' @return A list of Primo results with the following elements:
#' \tabular{ll}{
#' \code{post_prob} \tab matrix of posterior probabilities
#' (each column corresponds to an association pattern).\cr
#' \code{pis} \tab vector of estimated proportion of SNPs
#' belonging to each association pattern.\cr
#' \code{D_mat} \tab matrix of densities under each association pattern.\cr
#' \code{Gamma} \tab correlation matrix.\cr
#' }
#'
#' \itemize{
#' \item If the results were originally from the \eqn{t}-statistic version,
#' the list will additionally contain:
#' \tabular{ll}{
#' \code{Tstat_mod} \tab matrix of moderated t-statistics.\cr
#' \code{V_mat} \tab matrix of scaling factors under the alternative distribution.\cr
#' \code{mdf_sd_mat} \tab matrix of standard deviation adjustment according to
#'  moderated degrees of freedom: df/(df-2).\cr
#' \code{prior_df} \tab vector of the prior degrees of freedom for each marginal distribution.\cr
#' \code{prior_var} \tab vector of the prior variance estimators for each marginaldistribution.\cr
#' \code{unscaled_var} \tab vector of the unscaled variance priors on non-zero coefficients
#' for each marginal distribution.
#'  }
#'
#' \item If the results were originally from the \eqn{p}-value version,
#' the list will additionally contain:
#' \tabular{ll}{
#' \code{chi_mix} \tab matrix of \eqn{-2\log(P)}{-2*log(P)}-values.\cr
#' \code{A} \tab vector of scaling factors under the alternative distributions.\cr
#' \code{df_alt} \tab vector of degrees of freedom approximated for the alternative distributions.\cr
#'  }
#'  }
#'
#' @export
#'
subset_Primo_obj <- function(Primo_obj,idx){
  Primo_obj$post_prob <- Primo_obj$post_prob[idx,]
  Primo_obj$D_mat <- Primo_obj$D_mat[idx,]

  ## t-statistic version
  if(!is.null(Primo_obj$Tstat_mod)) Primo_obj$Tstat_mod <- Primo_obj$Tstat_mod[idx,]
  if(!is.null(Primo_obj$V_mat)) Primo_obj$V_mat <- Primo_obj$V_mat[idx,]
  if(!is.null(Primo_obj$mdf_sd_mat)) Primo_obj$mdf_sd_mat <- Primo_obj$mdf_sd_mat[idx,]

  ## p-value version
  if(!is.null(Primo_obj$chi_mix)) Primo_obj$chi_mix <- Primo_obj$chi_mix[idx,]

  return(Primo_obj)
}

#' Append a Primo object.
#'
#' Append two sets of Primo results. The function assumes that the two sets
#' share marginal distribution parameters and a common correlation structure.
#'
#' @param Primo_obj1 list of results returned by Primo (from the function
#' \code{\link{Primo_tstat}}, \code{\link{Primo_pval}}, or \code{\link{Primo_modT}}).
#' @param Primo_obj2 list of results returned by Primo.
#'
#' @inherit subset_Primo_obj return
#'
#' @export
#'
append_Primo_obj <- function(Primo_obj1,Primo_obj2){

  if(ncol(Primo_obj1$post_prob) != ncol(Primo_obj2$post_prob)) stop("Primo objects do not match (based on different numbers of studies).")

  Primo_obj1$post_prob <- rbind(Primo_obj1$post_prob,Primo_obj2$post_prob)
  Primo_obj1$D_mat <- rbind(Primo_obj1$D_mat,Primo_obj2$D_mat)

  ## t-statistic version
  if(!is.null(Primo_obj1$Tstat_mod)) Primo_obj1$Tstat_mod <- rbind(Primo_obj1$Tstat_mod,Primo_obj2$Tstat_mod)
  if(!is.null(Primo_obj1$V_mat)) Primo_obj1$V_mat <- rbind(Primo_obj1$V_mat,Primo_obj2$V_mat)
  if(!is.null(Primo_obj1$mdf_sd_mat)) Primo_obj1$mdf_sd_mat <- rbind(Primo_obj1$mdf_sd_mat,Primo_obj2$mdf_sd_mat)

  ## p-value version
  if(!is.null(Primo_obj1$chi_mix)) Primo_obj1$chi_mix <- rbind(Primo_obj1$chi_mix,Primo_obj2$chi_mix)

  return(Primo_obj1)
}


#' Find the lead SNP for each phenotype in each region.
#'
#' Determine the lead SNP for each phenotype in each region
#' based on summary statistics. A data.table will be returned
#' containing lead SNP information for each region.
#'
#' @param data data.table. Each row will be a SNP-phenotype combination
#' with statistics necessary to determine the lead SNP in each phenotype region.
#' @param snp_col character string of the column name of the SNP.
#' @param pheno_cols character vector of the column names of the phenotypes.
#' @param stat_cols character vector of the column names of statistics to be
#' used to determine lead SNPs.
#' @param data_type character string denoting type of statistics being used. Must be
#' either "pvalue" or "tstat".
#' @param suffices character vector denoting suffices to use for the names of the
#' lead SNP columns (optional). If \code{NULL}, consecutive integers will be assigned.
#'
#' @return A data.table containing information about the lead SNPs and
#' associated statistics. The columns will be \code{pheno_cols} followed by two columns
#' for each phenotype: the name of the lead SNP and the value of the statistic for
#' that lead SNP in that phenotype.
#'
#' @export
#'
find_leadsnps <- function(data,snp_col,pheno_cols,stat_cols,data_type="pvalue",suffices=NULL){

  base::requireNamespace("data.table")

  if(!data.table::is.data.table(data)){
    warning("Converting argument `data` to data.table format.")
    data <- data.table::data.table(data)
  }

  data.table::setkeyv(data,pheno_cols)

  if(is.null(suffices)) suffices <- 1:length(stat_cols)

  if(data_type=="pvalue"){

    leadsnps_region <- NULL

    for(i in 1:length(stat_cols)){

      ## find SNP with minimum p-value for the current phenotype in the region
      topSNP_currPheno <- data %>% dplyr::group_by(.dots=pheno_cols) %>% dplyr::slice(which.min(get(stat_cols[i])))
      topSNP_currPheno <- data.table::data.table(topSNP_currPheno,key=pheno_cols)
      topSNP_currPheno<- subset(topSNP_currPheno, select=c(pheno_cols,snp_col,stat_cols[i]))
      colnames(topSNP_currPheno)[ncol(topSNP_currPheno)-1] <- paste0("leadSNP_",suffices[i])

      ## merge results
      if(is.null(leadsnps_region)){
        leadsnps_region <- topSNP_currPheno
      } else{
        leadsnps_region <- merge(leadsnps_region,topSNP_currPheno)
      }

    }

  } else if(data_type=="tstat"){

    leadsnps_region <- NULL

    for(i in 1:length(stat_cols)){

      ## find SNP with maximum abs(t-statistic) for the current phenotype in the region
      topSNP_currPheno <- data %>% dplyr::group_by(.dots=pheno_cols) %>% dplyr::slice(which.max(abs(get(stat_cols[i]))))
      topSNP_currPheno <- data.table::data.table(topSNP_currPheno,key=pheno_cols)
      topSNP_currPheno<- subset(topSNP_currPheno, select=c(pheno_cols,snp_col,stat_cols[i]))
      colnames(topSNP_currPheno)[ncol(topSNP_currPheno)-1] <- paste0("leadSNP_",suffices[i])

      ## merge results
      if(is.null(leadsnps_region)){
        leadsnps_region <- topSNP_currPheno
      } else{
        leadsnps_region <- merge(leadsnps_region,topSNP_currPheno)
      }

    }

  } else{

    stop("data_type must be either 'pvalue' or 'tstat'.")

  }

  return(leadsnps_region)
}

#' Estimate posterior probabilities for observations missing from original Primo analysis.
#'
#' For each observation (e.g. SNP), estimates the posterior probability for each association pattern.
#' Uses parameters estimated by a previous run of \code{\link{Primo_tstat}} or \code{\link{Primo_modT}}
#' to estimate probabilities for SNPs missing in one or more studies.
#' Utilizes parallel computing, when available.
#'
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param dfs vector or matrix of degrees of freedom.
#' @param trait_idx integer vector of the columns corresponding to non-missing phenotypes/studies.
#' @param mafs vector or matrix of minor allele frequencies (MAFs).
#' If \code{NULL}, standard errors will not be adjusted for MAF.
#' @param pis matrix (one-row) of the estimated proportion of observations
#' belonging to each association pattern.
#' @param Gamma correlation matrix.
#' @param prior_df vector of the prior degrees of freedom for each marginal distribution.
#' @param prior_var vector of the prior variance estimators for each marginal distribution.
#' @param unscaled_var vector of the unscaled variance priors on non-zero coefficients
#' for each marginal distribution.
#' @param par_size numeric value; specifies the number of workers for
#' parallel computing (1 for sequential processing).
#'
#'
#' @inherit Primo_tstat return
#'
#' @export
#'
Primo_missdata_tstat <- function(betas,sds,dfs,trait_idx,mafs=NULL,pis,Gamma,prior_df,prior_var,unscaled_var,par_size=1){
  m <- nrow(betas)
  d <- ncol(betas)
  # orig_d <- log(length(pis),2)
  orig_d <- log(ncol(pis),2)

  miss_idx <- (1:orig_d)[-trait_idx]

  ## extract data for traits not missing data
  d0 <- prior_df[trait_idx]
  s02 <- prior_var[trait_idx]
  v0 <- unscaled_var[trait_idx]
  Gamma <- Gamma[trait_idx,trait_idx]

  orig_Q <- Primo::make_qmat(1:orig_d)
  if(length(miss_idx) == 1){
    keep_patterns <- which(orig_Q[,miss_idx]==0)
  } else{
    keep_patterns <- which(rowSums(orig_Q[,miss_idx])==0)
  }
  # pis <- pis[keep_patterns]
  pis <- pis[,keep_patterns]


  ## estimate marginal density functions in limma framework
  density_list <- lapply(1:d, function(j){
    if(is.matrix(mafs)) mafs <- mafs[,j]
    if(is.matrix(dfs)) {
      d1 <- dfs[,j]
    } else d1 <- rep(dfs[j],m)

    ## account for MAF in variance calculations
    if (is.null(mafs)){
      v1 = rep(1,m)
      sigma2 <- sds[,j]^2
    } else{
      v1 <- 1/(2*mafs*(1-mafs))
      sigma2 <- sds[,j]^2*(2*mafs*(1-mafs))
    }

    ## rescale t-statistic (see Smyth, 2004)
    sg_tilde <- sqrt((d0[j]*s02[j]+d1*sigma2)/(d0[j]+d1))
    moderate.t <- betas[,j]/(sg_tilde*sqrt(v1))

    ## estimate null and alternative densities
    df_mod=d0[j]+d1
    scaler=sqrt(1+v0[j]/v1)

    return(list(Tstat_mod = moderate.t, df_mod=df_mod, scaler=scaler, prior_df=d0[j], prior_var=s02[j], unscaled_var=v0[j]))
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
  Q<-Primo::make_qmat(1:d)
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


  ## obtain posterior probabilities
  PP <- sweep(log(D_mat),2,log(pis),"+")
  # subtract maximum, then e^(B)/rowSums(B)
  PP<-PP-matrixStats::rowMaxs(PP)
  PP<-exp(PP)
  # use rowSums directly since columns are recycled
  PP<-PP/rowSums(PP)

  return(list(post_prob=PP, pis=pis, D_mat=D_mat, Gamma=Gamma, Tstat_mod=Tstat_mod, V_mat = V, mdf_sd_mat = mdf_sd_mat,
              prior_df=prior_df,prior_var=prior_var,unscaled_var=unscaled_var))
}



#' Estimate posterior probabilities for observations missing from original Primo analysis.
#'
#' For each SNP, estimates the posterior probability for each association pattern.
#' Uses parameters estimated by a previous run of \code{\link{Primo_pval}} or \code{\link{Primo_chiMix}}
#' to estimate probabilities for SNPs missing in one or more studies.
#' \eqn{P}-values from non-missing studies are used as input.
#' Utilizes parallel computing, when available.
#'
#' @param pvals matrix of \eqn{P}-values from test statistics.
#' @param trait_idx integer vector of the columns corresponding to non-missing phenotypes/studies.
#' @param pis matrix (one-row) of the estimated proportion of observations
#' belonging to each association pattern
#' @param Gamma correlation matrix.
#' @param A vector of scaling factors under the alternative distributions.
#' @param df_alt vector of degrees of freedom approximated for the alternative distributions.
#' @param par_size numeric value; specifies the number of workers for
#' parallel computing (1 for sequential processing).
#'
#'
#' @inherit Primo_pval return
#'
#' @export
#'
Primo_missdata_pval <- function(pvals,trait_idx,pis,Gamma,A,df_alt,par_size=1){
  m <- nrow(pvals)
  d <- ncol(pvals)
  # orig_d <- log(length(pis),2)
  orig_d <- log(ncol(pis),2)

  miss_idx <- (1:orig_d)[-trait_idx]

  ## convert p-values to mixture of chi-squared statistics: -2log(p)
  chi_mix<-(-2)*log(pvals)

  ## extract data for traits not missing data
  A <- A[trait_idx]
  df_alt <- df_alt[trait_idx]
  Gamma <- Gamma[trait_idx,trait_idx]

  ## determine which patterns to keep
  orig_Q <- Primo::make_qmat(1:orig_d)
  if(length(miss_idx) == 1){
    keep_patterns <- which(orig_Q[,miss_idx]==0)
  } else{
    keep_patterns <- which(rowSums(orig_Q[,miss_idx])==0)
  }
  # pis <- pis[keep_patterns]
  pis <- pis[,keep_patterns]

  ## computation of D_mat (densities under each pattern)
  Q<-Primo::make_qmat(1:d)
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

  ## bind densities together into matrix
  D_mat <- do.call("cbind",D_mat_byk)
  rm(D_mat_byk)


  ## obtain posterior probabilities
  PP <- sweep(log(D_mat),2,log(pis),"+")
  # subtract maximum, then e^(B)/rowSums(B)
  PP<-PP-matrixStats::rowMaxs(PP)
  PP<-exp(PP)
  # use rowSums directly since columns are recycled
  PP<-PP/rowSums(PP)

  return(list(post_prob=PP, pis=pis, D_mat=D_mat, Gamma=Gamma, chi_mix=chi_mix, A=A, df_alt=df_alt))
}


#' Collapse posterior probabilities based on number of traits of association.
#'
#' Combine the posterior probabilities of association patterns
#' according to the number of traits with non-null associations in each pattern.
#' Provides the posterior probability of
#' being associated with "at least n" number of traits under analysis
#' by summing over association patterns with at least "n" studies
#' coming from the alternative distribution. Can also require non-null
#' association with one or more traits and combine over other traits.
#'
#' @param post_prob matrix of posterior probabilities.
#' @param req_idx (optional) scalar or integer vector of trait(s) where a non-null
#' association is required. If \code{NULL}, no traits will be required.
#' @param prefix character string denoting prefix of column names for the results matrix.
#'
#' @return A numeric matrix of combined posterior probabilities.
#'
#' @export
#'
collapse_pp_num <- function(post_prob,req_idx=NULL,prefix="pp_ge_"){

  ## number of traits
  d <- log(ncol(post_prob),2)
  ## pattern matrix
  Q <- Primo::make_qmat(1:d)

  ## columns for "at least X traits"
  if(is.null(req_idx)){
    cols_geX <- lapply(1:d, function(x) which(rowSums(Q) >= x))
    names(cols_geX) <- paste0(prefix,1:d)
  } else if(length(req_idx)==1){
    cols_geX <- lapply(0:(d-length(req_idx)), function(x) which(Q[,req_idx]==1 & rowSums(Q[,-req_idx]) >= x))
    names(cols_geX) <- paste0(prefix,0:(d-length(req_idx)))
  } else {
    cols_geX <- lapply(0:(d-length(req_idx)), function(x) which(rowSums(Q[,req_idx])==length(req_idx) & rowSums(Q[,-req_idx]) >= x))
    names(cols_geX) <- paste0(prefix,0:(d-length(req_idx)))
  }

  ## sum over columns denoting "at least X traits" to obtain combined posterior probability
  pp_geX <- lapply(cols_geX, function(idx){
    if(length(idx)==1){
      return(post_prob[,idx])
    } else{
      return(rowSums(post_prob[,idx]))
    }
  })

  pp_geX <- do.call("cbind",pp_geX)
  colnames(pp_geX) <- names(cols_geX)

  return(pp_geX)
}

#' Collapse posterior probabilities of association for each individual trait.
#'
#' Combine the posterior probabilities for each individual trait according to the
#' association patterns from which that trait comes from the alternative distribution.
#' Provides the posterior probability of being associated with "at least trait x".
#' Can also require non-null association with one or more traits and combine
#' over patterns where both "x" and that/those trait(s) are from the alternative distribution.
#'
#' @inheritParams collapse_pp_num
#'
#' @inherit collapse_pp_num return
#'
#' @export
#'
collapse_pp_trait <- function(post_prob,req_idx=NULL,prefix="pp_"){

  ## total number of traits
  d <- log(ncol(post_prob),2)
  ## pattern matrix
  Q <- Primo::make_qmat(1:d)

  ## columns for "at least X"
  if(is.null(req_idx)){
    cols_X <- lapply(1:d, function(x) which(Q[,x] == 1))
    names(cols_X) <- paste0(prefix,1:d)
  } else if(length(req_idx)==1){
    cols_X <- lapply((1:d)[-req_idx], function(x) which(Q[,req_idx]==1 & Q[,x] == 1))
    names(cols_X) <- paste0(prefix,(1:d)[-req_idx])
  } else {
    cols_X <- lapply((1:d)[-req_idx], function(x) which(rowSums(Q[,req_idx])==length(req_idx) & Q[,x] == 1))
    names(cols_X) <- paste0(prefix,(1:d)[-req_idx])
  }

  ## sum over columns denoting "at least X" to obtain combined posterior probability
  pp_X <- lapply(cols_X, function(idx){
    if(length(idx)==1){
      return(post_prob[,idx])
    } else{
      return(rowSums(post_prob[,idx]))
    }
  })

  pp_X <- do.call("cbind",pp_X)
  colnames(pp_X) <- names(cols_X)

  return(pp_X)
}

