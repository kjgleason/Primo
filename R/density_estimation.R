#' Estimate parameters under the null and alternative densities
#'
#' Estimate the parameters of the null and alternative
#' density functions using moderated \eqn{t}-statistics.
#'
#' @param betas vector of coefficient estimates.
#' @param sds vector of standard errors (for coefficient estimates).
#' @param df scalar or vector of first degrees of freedom of the F-distribution
#' (usually, \eqn{N-p-1}).
#' @param alt_prop proportion of test statistics from the alternative density.
#' @param mafs vector of minor allele frequencies (MAFs) (optional).
#' If \code{NULL}, no adjustment will be made for MAF in variance calculations.
#' @param N scalar or vector of number of subjects (optional).
#' Should be specified if \code{mafs!=NULL}.
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{Tstat_mod} \tab matrix of moderated t-statistics.\cr
#' \code{df_mod} \tab vector of the moderated degrees of freedom.\cr
#' \code{scaler} \tab vector of the scaling factors for the moderated t-statistics under the alternative.\cr
#' \code{prior_df} \tab scalar of the prior degrees of freedom.\cr
#' \code{prior_var} \tab scalar of the prior variance estimator.\cr
#' \code{unscaled_var} \tab scalar of the unscaled variance prior on non-zero coefficients.
#' }
#'
#' @details Following Smyth (2004), the function calculates moderated \eqn{t}-statistics using
#' conditional posterior means of the variance. The moderated \eqn{t}-statistics are used to estimate
#' the null density under a \eqn{t}-distribution and
#' the alternative density under a scaled \eqn{t}-distribution, each with moderated degrees of freedom.
#'
#' @export
#'
estimate_densities_modT <- function(betas, sds, df, alt_prop, mafs=NULL, N=NULL){
  m <- length(betas)

  ## if number of subjects is missing, add two to df (would be correct for simple linear regression; otherwise would be underestimate)
  if(is.null(N)) N <- df + 2

  ## account for MAF in variance calculations
  if (is.null(mafs)){
    v1 = rep(1,m)
    sigma2 <- sds^2
  } else{
    v1 <- 1/(2*mafs*(1-mafs)*N)
    sigma2 <- sds^2*(2*mafs*(1-mafs)*N)
  }

  ## estimate moments of scaled F-distribution using method of Smyth (2004)
  d1=df
  xx <- limma::fitFDist(sigma2,d1)
  s02 <- xx$scale; d0 <- xx$df2
  ## rescale t-statistic (see Smyth, 2004)
  sg_tilde <- sqrt((d0*s02+d1*sigma2)/(d0+d1))
  moderate.t <- betas/(sg_tilde*sqrt(v1))

  ## warn user about using low number of statistics to estimate alternative density
  pM <- alt_prop*length(betas)
  if(pM/2 < 30) warning("The specified proportion of alternative statistics yields a low total count.")

  ## estimate null and alternative densities
  df_mod=d0+d1
  v0 <- limma::tmixture.vector(moderate.t, sqrt(v1),df_mod,proportion=alt_prop,v0.lim=NULL)
  scaler=sqrt(1+v0/v1)

  return(list(Tstat_mod = moderate.t, df_mod=df_mod, scaler=scaler, prior_df=d0, prior_var=s02, unscaled_var=v0))
}

#' Estimate parameters under the alternative density
#'
#' Estimate the parameters of the alternative
#' density function using a mixture of chi-squared distributions of
#' \eqn{-2\log(P)}{-2*log(P)}-values.
#'
#' @param pvals vector of \eqn{P}-values.
#' @param alt_prop proportion of \eqn{P}-values used in estimating alternative densities.
#' @param method_moments logical, denoting whether to estimate scale and degree of freedom
#' parameters using Method of Moments (\code{TRUE}) or optimization (\code{FALSE}).
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{chi_mix} \tab vector of \eqn{-2\log(P)}{-2*log(P)}-values.\cr
#' \code{scaler} \tab estimated scaling factor for the alternative distribution.\cr
#' \code{df_alt} \tab estimated degrees of freedom for the alternative distribution.\cr
#' }
#'
#' @details The function estimates densities under the null and under the alternative
#' given a vector of \eqn{P}-values from test statistics and the marginal probability of
#' coming from the alternative distribution. Under the null hypothesis,
#' \code{-2*log()} transformed \eqn{P}-values follow a chi-squared distribution with 2 degrees of freedom.
#' Under the alternative, the function assumes that the \code{-2*log()} transformed \eqn{P}-values
#' follow a scaled chi-squared distribution with unknown scale parameter and
#' unknown degrees of freedom.
#'
#' When \code{method_moments=TRUE}, the two unknown parameters are estimated by using the first and
#' second moments of the transformed \eqn{P}-values given the known proportion of alternatives
#' (specified in \code{alt_prop}). The function solves a theoretical formula of the first and
#' second moments as functions of the scale parameter and degrees of freedom. Note that
#' the alternative density estimation may fail if \code{alt_prop} is small and
#' \code{method_moments=TRUE}.
#'
#' When \code{method_moments=FALSE}, the two unknown parameters are estimated by minimizing the differences
#' between the \eqn{P}-values given the parameters and the nominal \eqn{P}-values based on each statistic's rank,
#' using optimization algorithms.
#'
#' @export
#'
estimate_densities_pval <- function(pvals, alt_prop, method_moments=FALSE){

  ## transform to chi-squared mixture (df=2 under null; scaled with unknown df under alternative)
  chi_mix<-(-2)*log(pvals)

  if(method_moments){
    if(alt_prop < 0.01) warning("magnitude of alt_prop is small; alternative density estimation may fail")

    ## Use the first and second moments to solve the scale parameter and the degrees of freedom under the alternative
    prod<-(mean(chi_mix)-(1-alt_prop)*2)/alt_prop
    ## store scale parameter in a_alt
    a_alt<-(mean((chi_mix-mean(chi_mix))^2)-(1-alt_prop)*2*2-(1-alt_prop)*alt_prop*(prod-2)^2)/(alt_prop*2*prod)
    ## store degrees of freedom in df_alt
    df_alt<-prod/a_alt

  } else{

    optim_dat <- list(chi_mix=sort(chi_mix, decreasing=T), alt_prop=alt_prop)
    init1 <- 2
    init2 <- 3

    ## run global optimization to identify correct "neighborhood" of optimum
    global_res <- nloptr::nloptr(x0=c(init1,init2),eval_f=chiMix_pDiff,lb=c(1,2),ub=c(100,100),
                                 opts=list(algorithm="NLOPT_GN_DIRECT",maxeval=500),
                                 data=optim_dat, sorted=T)
    init1 <- global_res$solution[1]
    init2 <- global_res$solution[2]

    ## refine optimum using a local optimization algorithm
    optim_res <- nloptr::nloptr(x0=c(init1,init2),eval_f=chiMix_pDiff,lb=c(1,2),
                                opts=list(algorithm="NLOPT_LN_COBYLA",maxeval=500),
                                data=optim_dat, sorted=T)

    ## store scale parameter in a_alt
    a_alt<-optim_res$solution[1]
    ## store degrees of freedom in df_alt
    df_alt<-optim_res$solution[2]
  }

  return(list(chi_mix=chi_mix, scaler=a_alt, df_alt=df_alt))
}

#' Difference from nominal \eqn{P}-values of a chi-squared mixture.
#'
#' Match observed \eqn{P}-values to nominal \eqn{P}-values by rank and calculate total difference. Objective
#' function of a chi-squared mixture is used. This function is optimized to estimate the unknown
#' scale and degree of freedom parameters for the scaled chi-squared distribution under the
#' alternative hypothesis.
#'
#' @param par vector of two parameters for the alternative distribution:
#' scale parameter [1] and degrees of freedom [2].
#' @param data list with two elements: \code{chi_mix} (vector) and \code{alt_prop} (scalar).
#' @param sorted logical, denoting whether \code{data$chi_mix} is sorted in decreasing order.
#'
#' @return A numeric value of the total (absolute) difference between observed p-values
#' (given the parameters) and nominal p-values.
#'
#' @details The argument \code{data} should be a two-element list. The first element, named \code{chi_mix}, is
#' a vector of the observed values from the chi-squared mixture. The second element, named \code{alt_prop},
#' is the proportion of statistics that come from the alternative distribution.
#'
#' The order of parameters is important in \code{par}. The scale parameter should be in the first position;
#' the degrees of freedom should be in the second position.
#'
#' @export
#'
chiMix_pDiff <- function(par, data, sorted=F){
  chi_mix <- data$chi_mix
  alt_prop <- data$alt_prop

  ## total number of statistics
  M <- length(chi_mix)

  ## number of statistics to use in estimation
  pM <- alt_prop * M

  ## warn user about using low number of statistics to estimate alternative density
  if(pM/2 < 30){
    if(pM/2 >= 20){
      warning("The specified proportion of alternative statistics yields a low total count.")
    } else warning(paste("The specified proportion of alternative statistics yields a low total count.",
                         "\nFitting density to top 20 test-statistics."))
  }

  ##===== Fit top (p*M/2) statistics to objective function (mixture distribution) =====##
  ## decreasing ranks of chi_mix (restricted to those being used to estimate alternative density parameters)
  r <- 1:max(pM/2,20)

  if(!sorted) chi_mix <- sort(chi_mix, decreasing=T)
  chi_mix <- chi_mix[r]

  ## par[1] holds A, the scale factor; par[2] holds d, the degrees of freedom
  obj_sum <- sum(abs(alt_prop*stats::pgamma(chi_mix,shape=par[2]/2,scale=2*par[1],lower.tail=F) + (1-alt_prop)*stats::pchisq(chi_mix,2,lower.tail=F) - (r-0.5)/M))

  return(obj_sum)
}

#' Estimate parameters under the null and alternative densities
#'
#' Estimate the parameters of the null and alternative
#' density functions.
#' This function calls either \code{\link{estimate_densities_modT}} or
#' \code{\link{estimate_densities_pval}} for density estimation using
#' moderated \eqn{t}-statistics or \eqn{P}-values, respectively.
#'
#' @param betas vector of coefficient estimates.
#' @param sds vector of standard errors (for coefficient estimates).
#' @param df scalar or vector of first degrees of freedom of the F-distribution.
#' @param pvals vector of \eqn{P}-values.
#' @param alt_prop proportion of test statistics from the alternative density.
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param N scalar or vector of number of subjects.
#' @param use_tstats logical indicating whether to use moderated \eqn{t}-statistics (\code{TRUE})
#' or \eqn{P}-values (\code{FALSE}).
#'
#' @return See \code{\link{estimate_densities_modT}} or \code{\link{estimate_densities_pval}}.
#'
#' @export
#'
#'
estimate_densities <- function(betas=NULL, sds=NULL, df=NULL, pvals=NULL, alt_prop, mafs=NULL, N=NULL, use_tstats=TRUE){

  if(use_tstats){
    if(is.null(betas) | is.null(sds) | is.null(df)){
      if(!is.null(pvals)) {
          warning("use_tstats=TRUE requires non-null values for betas, sds, mafs and df; using p-values for density estimation.")
          myDens <- estimate_densities_pval(pvals=pvals, alt_prop=alt_prop)
      } else stop("use_tstats=TRUE requires non-null values for betas, sds, and df.")
    } else myDens <- estimate_densities_modT(betas=betas, sds=sds, df=df, alt_prop=alt_prop, mafs=mafs, N=N)
  } else{
    myDens <- estimate_densities_pval(pvals=pvals, alt_prop=alt_prop)
  }
  return(myDens)
}

#' Create association pattern matrix
#'
#' Make matrix denoting association patterns (\eqn{Q}) that indicates null/alternative densities
#' in the association patterns.
#'
#' @param grp vector of group identifiers (common to use \code{seq(1,d)}, for \eqn{d} groups).
#' @param name_cols logical indicating whether to name columns using group identifiers.
#'
#' @return A matrix with all binary combinations of the groups. The value in row \eqn{k}
#' and column \eqn{j} denotes whether statistics from study \eqn{j} come from
#' the null (\eqn{0}) or alternative (\eqn{1}) distribution
#' under the \eqn{k}-th association pattern.
#'
#' @export
#'
make_qmat <- function(grp, name_cols=FALSE){

  ng <- length(grp)
  ngi <- rep(1,ng)
  Q <- rep(0,ng)

  for (qi in 1:ng){
    oo <- utils::combn(1:ng,qi)
    Q <- cbind(Q, apply(oo,2,function(x,ng,ngi) {
      vv <- rep(0,ng)
      vv[x] <- 1
      oo <- rep(vv,ngi)
      return(oo)
    }, ng=ng, ngi=ngi))
  }
  Q <- t(Q)

  dimnames(Q) <- NULL

  if(name_cols) colnames(Q) <- grp

  return(Q)
}
