#' Estimate Densities Under the Null and Alternative Densities
#'
#' For each observation, estimate the density under the null and under the alternative
#' hypotheses using moderated t-statistics.
#'
#' @param betas vector of coefficient estimates.
#' @param sds vector of standard errors (for coefficient estimates).
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param df first degrees of freedom of the F-distribution (usually, the number of subjects/observations).
#' @param alt_proportion proportion of test-statistics used in estimating alternative densities.
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{Tstat_m} \tab matrix of moderated t-statistics\cr
#' \code{D0} \tab matrix of densities calculated under the null distribution\cr
#' \code{D1} \tab matrix of densities calculated under the alternative distribution\cr
#' }
#'
#' @details Following Smyth (2004), the function calculates moderated t-statistics using
#' conditional posterior means of the variance. The moderated t-statistics are used to estimate
#' the null density under a t-distribution with estimated degrees of freedom and
#' the alternative density under a scaled t-distribution.
#'
#' @export
#'
estimate_densities_modT <- function(betas, sds, mafs, df, alt_proportion){

  # account for MAF in variance calculations
  vg = 1/(2*mafs*(1-mafs))
  sigma2 <- sds^2*(2*mafs*(1-mafs))

  # estimate moments of scaled F-distribution using method of Smyth (2004)
  d1=df
  xx <- limma::fitFDist(sigma2,d1)
  s02 <- xx$scale; n0 <- xx$df2
  # rescale t-statistic (see Smyth, 2004)
  sg_tilde <- sqrt((n0*s02+d1*sigma2)/(n0+d1))
  moderate.t <- betas/(sg_tilde*sqrt(vg))

  # defunct: old way of estimating proportion of null
  #pt <- 2*(1-pt(abs(moderate.t), df=d1+n0))
  #pi0 <- 1
  #try({pi0 <- qvalue(pt)$pi0}, silent=TRUE)

  # estimate null and alternative densities
  D0 <- dt(moderate.t, df=d1+n0)
  v0 <- limma::tmixture.vector(moderate.t, sqrt(vg),d1+n0,proportion=alt_proportion,v0.lim=NULL)
  scaler=sqrt(1+v0/vg)
  D1 <- metRology::dt.scaled(moderate.t,df=d1+n0,mean=0,sd=scaler)

  return(list(Tstat_m = moderate.t, D0=D0, D1=D1))
}

#' Estimate Densities Under the Null and Alternative Densities
#'
#' For each observation, estimate the density under the null and under the alternative
#' hypotheses using p-values.
#'
#' @param pvals vector of coefficient estimates.
#' @param alt_proportion proportion of pvalues used in estimating alternative densities.
#' Note that the alternative density estimation may fail
#' if the value is too small.
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{D0} \tab matrix of densities calculated under the null distribution\cr
#' \code{D1} \tab matrix of densities calculated under the alternative distribution\cr
#' }
#'
#' @details The function estimates densities under the null and under the alternative
#' given a vector of p-values from test statistics and the marginal probability of
#' coming from the alternative distribution. Under the null hypothesis,
#' \code{-2*log()} transformed p-values follow a chi-square distribution with 2 degrees of freedom.
#' Under the alternative, the function assumes that the \code{-2*log()} transformed p-values
#' follow a scaled chi-square distribution with unknown scale parameter and
#' unknown degrees of freedom. Those two parameters are solved by using the first and
#' second moments of the transformed p-values and the known proportion of alternatives
#' (specified in \code{alt_proportion}). The function solves a theoretical formula of the first and
#' second moments as functions of the scale parameter and degrees of freedom.
#'
#' @export
#'
estimate_densities_pval <- function(pvals, alt_proportion){

  if(alt_proportion < 0.05) warning("magnitude of alt_proportion is small; alternative density estimation may fail")

  ##Transform to chi square; under the null, the transformed p-values follow a chi square distribution with df 2.
  chi_mix<-(-2)*log(pvals)

  ##Using the first and second moments to solve the scale parameter and the degrees of freedom under the alternative
  prod<-(mean(chi_mix)-(1-alt_proportion)*2)/alt_proportion
  ## store scale parameter in a_alt
  a_alt<-(mean((chi_mix-mean(chi_mix))^2)-(1-alt_proportion)*2*2-(1-alt_proportion)*alt_proportion*(prod-2)^2)/(alt_proportion*2*prod)
  ## store degrees of freedom in df_alt
  df_alt<-prod/a_alt

  ##Density under the null
  D0<-dchisq(chi_mix,df=2)
  ##Density under the alternative
  D1<-dchisq(chi_mix/a_alt,df=df_alt)/a_alt

  return(list(D0=D0, D1=D1))
}

#' Estimate Densities Under the Null and Alternative Densities
#'
#' For each observation, estimate the density under the null and under the alternative hypotheses.
#' This function is a wrapper to call either \code{estimate_densities_modT()}or
#' \code{estimate_densities_pval()} for density estimation using moderated t-statistics or p-values, respectively.
#'
#' @param pvals vector of p-values.
#' @param betas vector of coefficient estimates.
#' @param sds vector of standard errors (for coefficient estimates).
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param df first degrees of freedom of the F-distribution.
#' @param alt_proportion proportion of test-statistics used in estimating alternative densities.
#' @param use_tstats logical; when true densities are calculated using moderated t-statistics
#' (use p-values when false)
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{Tstat_m} \tab matrix of moderated t-statistics (set to NULL for p-value method)\cr
#' \code{D0} \tab matrix of densities calculated under the null distribution\cr
#' \code{D1} \tab matrix of densities calculated under the alternative distribution\cr
#' }
#'
#' @export
#'
#'
estimate_densities <- function(pvals=NULL, betas=NULL, sds=NULL, mafs=NULL, df=NULL, alt_proportion, use_tstats=TRUE){

  if(use_tstats){
    if(is.null(betas) | is.null(sds) | is.null(mafs) | is.null(df)){
      if(!is.null(pvals)) {
          warning("use_tstats=TRUE requires non-null values for betas, sds, mafs and df; using p-values for density estimation.")
          myDens <- estimate_densities_pval(pvals, alt_proportion)
      } else stop("use_tstats=TRUE requires non-null values for betas, sds, mafs and df.")
    } else myDens <- estimate_densities_modT(betas, sds, mafs, df, alt_proportion)
  } else{
    myDens <- estimate_densities_pval(pvals, alt_proportion)
    myDens <- list(Tstat_m=NULL,D0=myDens$D0,D1=myDens$D1)
  }
  return(myDens)
}



#' Estimate Posterior Probabilities of Configurations
#'
#' For each SNP, estimates the posterior probability for each configuration.
#' Utilizes parallel computing, when available.
#'
#' @param pvals matrix of p-values from test statistics.
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param dfs vector of degrees of freedom.
#' @param alt_proportions vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param tol numerical value; specifies tolerance threshold for convergence.
#' @param par_size numerical value; specifies the number of CPUs/cores/processors for
#' parallel computing (0 for sequential processing).
#' @param density_list (optional) list of densities estimated by \code{estimate_densities()}.
#' If provided by user, must follow the format specified in Details.
#'
#' @return A list with the following elements:
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
#' @details Either \code{pvals} are all four of \code{betas}, \code{sds}, \code{mafs} and \code{dfs}
#' must be specified (i.e. not \code{NULL}).
#'
#' The following are additional details describing the input arguments:
#' \tabular{ll}{
#' \code{pvals} \tab  p-values should be an LxD matrix: L = number of SNPs ("Loci");
#' D = number data sources.\cr
#' \code{betas} \tab  beta coefficients should be an LxD matrix: L = number of SNPs ("Loci");
#' D = number data sources.\cr
#' \code{sds} \tab standard errors should be an LxD matrix: L = number of SNPs ("Loci");
#' D = number data sources.\cr
#' \code{mafs} \tab minor allele frequencies (MAF) should be vector with length matching
#' the number of SNPs/Loci\cr
#'  \tab (same order as rows in \code{betas} and \code{sds}).\cr
#' \code{dfs} \tab vector containing number of subjects in each dataset\cr
#'  \tab (same order as columns in \code{betas} and \code{sds})\cr
#' \code{alt_proportions} \tab proportion of test-statistics used to estimate
#' the alternative distribution for each data source\cr
#'  \tab (same order as columns in \code{betas} and \code{sds})\cr
#' \code{density_list} \tab list must be 3 elements, named \code{Tstat_m},\code{D0}, and \code{D1}
#'  -- outputs from \code{estimate_densities()}\cr
#'  \tab (each element is a matrix, with row matching observations of \code{betas} and \code{sds},\cr
#'  \tab and columns in the same order as columns in \code{betas} and \code{sds})
#' }
#'
#' @export
#'
estimate_config <- function(pvals=NULL, betas=NULL, sds=NULL, mafs=NULL, dfs=NULL, alt_proportions, tol=1e-3, par_size=0, density_list=NULL){

  # store dimensions of test statistics
  if(!is.null(betas)){
    m <- nrow(betas)
    d <- ncol(betas)
  } else{
    m <- nrow(pvals)
    d <- ncol(pvals)
  }

  ## estimate null and alternate densities for each column/study
  if(is.null(density_list)){
    Tstat_m <- D0 <- D1 <- NULL
    if(!is.null(betas)){
      for (j in 1:d){
        temp <- estimate_densities(betas=betas[,j],sds=sds[,j],mafs=mafs,df=dfs[j],alt_proportion = alt_proportions[j])
        Tstat_m <- cbind(Tstat_m, temp$Tstat_m)
        D0 <- cbind(D0, temp$D0)
        D1 <- cbind(D1, temp$D1)
      }
    } else{
      for (j in 1:d){
        temp <- estimate_densities(pvals=pvals[,j],alt_proportion = alt_proportions[j],use_tstats=F)
        Tstat_m <- cbind(Tstat_m, temp$Tstat_m)
        D0 <- cbind(D0, temp$D0)
        D1 <- cbind(D1, temp$D1)
      }
    }
  }else{
    Tstat_m <- density_list$Tstat_m
    D0 <- density_list$D0
    D1 <- density_list$D1
  }

  if (par_size>0){
    cat("\nCreating cluster...\n")
    cl <- parallel::makeCluster(par_size)
  } else cl <- NULL

  # fix configuration matrix and number of patterns
  Q <- make_qmat(1:d)
  n_pattern <- nrow(Q)
  # initialize values for proportions coming from each configuration
  curpi<- c(0.80, rep((1-0.80)/(n_pattern-1),n_pattern-1))
  diff<-1
  numiters<-1
  itermat<-curpi

  cat("\nIterating...\n")

  # EM algorithm
  # maximum vector size in R is 2^31-1 (need to process in chunks if # elements in matrix exceeds)
  if(m*as.double(n_pattern) <= 2^31-1){
    while(diff>tol){
      numiters<-numiters+1
      curb<-e_step(curpi, Q=Q, D_0=D0, D_1=D1)
      # temporary solution to precision problems caused by |t-statistics| > 30
      curb[which(is.na(curb))] <- 0
      curpi<-m_step(curb)
      itermat<-rbind(itermat,curpi)
      diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
      if (!(numiters %% 10)) cat("\nIteration:",numiters,"; Change:",diff,"\nPi-hat:",curpi,"\n")
    }
  } else{
    while(diff>tol){
      numiters<-numiters+1
      # process large matrix in chunks
      num_chunks <- ceiling((m*as.double(n_pattern))/(2^31-1))
      Drow_chunks <- split(1:m, ceiling(seq_along(1:m)/(m/num_chunks)))
      # e-step, in chunks
      curb_colsums <- 0
      for(ch in 1:length(Drow_chunks)){
        D_rows <-  Drow_chunks[[ch]]
        curb_chunk <- e_step(curpi,Q,D0[D_rows,],D1[D_rows,])
        # temporary solution to precision problems caused by |t-statistics| > 30
        curb_chunk[which(is.na(curb_chunk))] <- 0
        curb_colsums <- curb_colsums + colSums(curb_chunk)
      }
      curpi<-matrix((curb_colsums+1)/(m+n_pattern),nrow=1)
      itermat<-rbind(itermat,curpi)
      diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
      if (!(numiters %% 10)) { cat("\nIteration:",numiters,"; Change:",diff,"\nPi-hat:",curpi,"\n")
      } else cat("\nIteration:",numiters,"; Change:",diff,"\n")
    }
  }


  if (par_size>0) parallel::stopCluster(cl)

  cat("\nIteration:",numiters,"; Change:",diff,"\nPi-hat:",curpi,"\n\n")

  if(m*as.double(n_pattern) <= 2^31-1){
    return(list(post_prob = curb, config_prop = curpi, Tstat_m = Tstat_m, D0=D0, D1=D1,
                alt_proportions=alt_proportions, tol=tol))
  } else{
    curb <- NULL
    exit_m_step <-function(curpi,Q,D0,D1){
      t_Q <- t(Q)
      # consider using `bigmemory` package for following matrix
      curb <- log(D0) %*% (1-t_Q) + log(D1) %*% t_Q
      curb <- sweep(curb,2,log(curpi),"+")
      # subtract minimum, then e^(B)/rowSums(B)
      curb<-curb-matrixStats::rowMins(curb)
      curb<-exp(curb)
      # can use rowSums directly since columns are recycled
      curb<-curb/matrixStats::rowSums2(curb)
      return(curb)
    }
    try(curb<-exit_m_step(curpi,Q,D0,D1))
    return(list(post_prob = curb, config_prop = curpi, Tstat_m = Tstat_m, D0=D0, D1=D1,
                alt_proportions=alt_proportions, tol=tol))
  }
}


#' Make Configuration Matrix
#'
#' Make configuration matrix (Q) that indicates null/alternative density patterns.
#'
#' @param grp vector of group identifiers (common to use seq(1,N), for N groups)
#' @param name.cols logical; whether or not to name columns using group identifiers
#'
#' @return \code{make_qmat} returns a matrix with all binary combinations of the groups
#'
#' @export
#'
make_qmat <- function(grp, name.cols=F){
  ## expand.grid simpler, but order of groups not as logical
  # ng <- length(grp)
  # Q <- as.matrix(expand.grid(replicate(ng,0:1,simplify=F)))
  # if(name.cols) colnames(Q) <- grp
  # else colnames(Q) <- NULL

  ng <- length(grp)
  ngi <- rep(1,ng)
  Q <- rep(0,ng)

  for (qi in 1:ng){
    oo <- combn(1:ng,qi)
    Q <- cbind(Q, apply(oo,2,function(x,ng,ngi) {
      vv <- rep(0,ng)
      vv[x] <- 1
      oo <- rep(vv,ngi)
      return(oo)
    }, ng=ng, ngi=ngi))
  }
  Q <- t(Q)

  if(name.cols) colnames(Q) <- grp

  return(Q)
}
