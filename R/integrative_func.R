#' Estimate Posterior Probabilities of Configurations
#'
#' For each SNP, estimates the posterior probability for each configuration.
#' Utilizes parallel computing, when available.
#'
#' @param betas matrix of coefficient estimates.
#' @param sds matrix of standard errors (for coefficient estimates).
#' @param mafs vector of minor allele frequencies (MAFs).
#' @param dfs vector of degrees of freedom.
#' @param alt_proportions vector of the proportions of test-statistics used in estimating
#' alternative densities.
#' @param tol numerical value; specifies tolerance threshold for convergence.
#' @param par_size numerical value; specifies the number of CPUs/cores/processors for
#' parallel computing
#' (0 for sequential processing).
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
#' @details The following are additional details describing the input arguments:
#' \tabular{ll}{
#' \code{betas} \tab  coefficients should be an LxD matrix: L = number of SNPs ("Loci");
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
#'  \tab (same order as columns in \code{betas} and \code{sds})
#' }
#'
#' @export
#'
estimate_config <- function(betas, sds, mafs, dfs, alt_proportions, tol=1e-3, par_size=0){

  # store dimensions of test statistics
  m <- nrow(betas)
  d <- ncol(betas)
  # account for MAF in variance calculations
  vg = 1/(2*mafs*(1-mafs))
  sigma2 <- sds^2*matrix(rep(2*mafs*(1-mafs), each=d), byrow=T, ncol=d)

  Tstat_m <- D0 <- D1 <- NULL

  for (j in 1:d){
    # estimate moments of scaled F-distribution using method of Smyth (2004)
    d1=rep(dfs[j],m)
    xx <- limma::fitFDist(sigma2[,j],d1)
    s02 <- xx$scale; n0 <- xx$df2
    # rescale t-statistic (see Smyth, 2004)
    sg_tilde <- sqrt((n0*s02+d1*sigma2[,j])/(n0+d1))
    moderate.t <- betas[,j]/(sg_tilde*sqrt(vg))
    Tstat_m <- cbind(Tstat_m, moderate.t)

    # defunct: old way of estimating proportion of nulls
    #pt <- 2*(1-pt(abs(moderate.t), df=d1+n0))
    #pi0 <- 1
    #try({pi0 <- qvalue(pt)$pi0}, silent=TRUE)

    # estimate null and alternative densities
    D0 <- cbind(D0, dt(moderate.t, df=d1+n0))
    v0 <- limma::tmixture.vector(moderate.t, sqrt(vg),d1+n0,proportion=alt_proportions[j],v0.lim=NULL)
    scaler=sqrt(1+v0/vg)
    D1 <- cbind(D1, metRology::dt.scaled(moderate.t,df=d1+n0,mean=0,sd=scaler))
  }

  if (par_size>0){
    cat("\nCreating cluster...\n")
    cl <- parallel::makeCluster(par_size)
  } else cl <- NULL

  # fix configuration matrix and number of patterns
  Q <- makeQ(1:d)
  n_pattern <- nrow(Q)
  # initialize values for proportions coming from each configuration
  curpi<- c(0.80, rep((1-0.80)/(n_pattern-1),n_pattern-1))
  diff<-1
  numiters<-1
  itermat<-curpi

  cat("\nIterating...\n")

  # EM algorithm
  while(diff>tol){
    numiters<-numiters+1
    curb<-e_step(curpi, n_obs=m, n_pattern=n_pattern, q=Q, density_0=D0, density_1=D1,cl=cl)
    curpi<-m_step(curb, n_obs=m, n_pattern=n_pattern)
    itermat<-rbind(itermat,curpi)
    diff<-sum(abs(itermat[numiters,]-itermat[numiters-1,]))/sum(itermat[numiters-1,])
    if (!(numiters %% 10)) cat("\nIteration:",numiters,"Diff:",diff,"\nPi-hat:",curpi,"\n")
  }

  if (par_size>0) parallel::stopCluster(cl)

  cat("\nIteration:",numiters,"; Change:",diff,"\nPi-hat:",curpi,"\n\n")

  return(list(post_prob = curb, config_prop = curpi, Tstat_m = Tstat_m, D0=D0, D1=D1))
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
