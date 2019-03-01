#' Perform fine-mapping on a specified variant.
#'
#' For a specified variant, re-estimate the posterior probabilities of association
#' patterns, conditioning on other specified variants. Returns the association
#' pattern with the highest posterior probability after fine-mapping.
#'
#' @param idx.snp scalar of the index of the genetic variant (e.g. row of \code{Tstat_mod})
#' on which one wants to perform fine-mapping
#' @param idx.leadsnps vector of indices of the leading snps (e.g. rows of \code{Tstat_mod})
#' on which one wants to condition during fine-mapping
#' @param LD_mat matrix of LD coefficients (\eqn{r^{2}}{r^2}).
#' Rows and columns should match the order of (\code{idx.snp}, \code{idx.leadsnps}).
#' @param Primo_obj list returned by running the \eqn{t}-statistic version
#' of Primo (i.e. \code{\link{Primo_tstat}})
#'
#' @return The numerical value corresponding to the association pattern
#' with the highest posterior probability following fine-mapping adjustment.
#' The value returned, \eqn{k}, corresponds to the \eqn{k}-th column of \code{pis}
#' from the Primo output, and the \eqn{k}-th row of the \eqn{Q} matrix produced
#' by \code{\link{make_qmat}}.
#'
#'
#' @export
#'
fine_map<-function(idx.snp,idx.leadsnps,LD_mat,Primo_obj){
  n_leadsnps<-length(idx.leadsnps)
  zi<-Primo_obj$Tstat_mod[c(idx.snp,idx.leadsnps),]
  d<-ncol(zi)
  Q <- make_qmat(1:d)
  V <- Primo_obj$V_mat[c(idx.snp,idx.leadsnps),]
  mdf_sd<-Primo_obj$mdf_sd_mat[c(idx.snp,idx.leadsnps),]
  post_prob<-Primo_obj$post_prob
  pis<-Primo_obj$pis
  v2<-NULL
  Gamma<-Primo_obj$Gamma
  for(i in 1:n_leadsnps){
    q2<-Q[which.max(post_prob[idx.leadsnps[i],]),]
    v2<-rbind(v2,(V[(i+1),]%*%diag(q2)+ matrix(1, nrow=1,ncol=d)%*%diag( 1-q2))*matrix(mdf_sd[(i+1),],ncol=d))
  }

  CD_mat=matrix(NA, nrow=1, ncol=2^d)

  for(k in 1:2^d){
    v_k <-  (V[1,]%*%diag(Q[k,]) + matrix(1, nrow=1,ncol=d)%*%diag( 1-Q[k,]))*matrix(mdf_sd[1,],ncol=d)
    cmean<-NULL
    csigma<-NULL
    for(j in 1:d){
      Var<-(c(v_k[j],v2[,j])%*%t(c(v_k[j],v2[,j])))*sqrt(LD_mat)
      cmean<-c(cmean,Var[1,(2:(1+n_leadsnps))]%*%solve(Var[(2:(1+n_leadsnps)),(2:(1+n_leadsnps))])%*%zi[(2:(1+n_leadsnps)),j])
      csigma<-c(csigma,Var[1,1] - Var[1,(2:(1+n_leadsnps))]%*%solve(Var[(2:(1+n_leadsnps)),(2:(1+n_leadsnps))])%*%as.matrix(Var[1,(2:(1+n_leadsnps))]))
    }
    sigma<-(sqrt(csigma)%*%t(sqrt(csigma)))*Gamma
    CD_mat[,k]<-mvtnorm::dmvnorm(zi[1,],mean = cmean,sigma =sigma,log=TRUE)
  }

  sp<-which.max(CD_mat+log(pis))

  # return(list(sp = sp))
  return(sp)
}
