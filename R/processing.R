#' Subset a Primo object.
#'
#' Subset results from Primo output based on a vector of indices.
#'
#' @param Primo_obj A list of results returned by Primo (from the function
#' \code{\link{Primo_tstat}}, \code{\link{Primo_pval}}, or \code{\link{Primo_ModT}}).
#' @param idx Integer vector of the indices to which to subset Primo results.
#'
#' @return List of Primo results with the following elements:
#' \tabular{ll}{
#' \code{post_prob} \tab matrix of posterior probabilities
#' (rows are SNPs; columns are association patterns)\cr
#' \code{pis} \tab vector of estimated proportion of SNPs
#' belonging to each association pattern\cr
#' \code{D_mat} \tab matrix of densities under each association pattern\cr
#' \code{Gamma} \tab correlation matrix\cr
#' }
#'
#' If the results were originally from the \eqn{t}-statistic version,
#' the list will additionally contain:
#' \tabular{ll}{
#' \code{Tstat_mod} \tab matrix of moderated t-statistics\cr
#' \code{V_mat} \tab matrix of scaling factors under the alternative distribution\cr
#' \code{mdf_sd_mat} \tab matrix of standard deviation adjustment according to
#'  moderated degrees of freedom: df/(df-2)\cr
#'  }
#'
#' If the results were originally from the \eqn{p}-value version,
#' the list will additionally contain:
#' \tabular{ll}{
#' \code{chi_mix} \tab matrix of \eqn{-2}log(\eqn{P})-values \cr
#' \code{A} \tab matrix of scaling factors under the alternative distribution\cr
#' \code{df_alt} \tab matrix of standard deviation adjustment according to
#'  moderated degrees of freedom: df/(df-2)\cr
#'  }
#'
#' @export
#'
subset_Primo_obj <- function(Primo_obj,idx){
  Primo_obj$Tstat_mod <- Primo_obj$Tstat_mod[idx,]
  Primo_obj$post_prob <- Primo_obj$post_prob[idx,]
  Primo_obj$D_mat <- Primo_obj$D_mat[idx,]
  Primo_obj$V_mat <- Primo_obj$V_mat[idx,]
  Primo_obj$mdf_sd_mat <- Primo_obj$mdf_sd_mat[idx,]

  return(Primo_obj)
}
