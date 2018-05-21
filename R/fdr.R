#' Calculate empirical False Discovery Rate (eFDR)
#'
#' Calculate the empirical False Discovery Rate (eFDR) given
#' observed results and results from permuted datasets. Takes
#' a vector of potential posterior probability thresholds and
#' calculates the eFDR for each threshold.
#'
#' @param true_res a list; results returned by \code{\link{estimate_config}} run
#' with real data
#' @param perm_res a list of lists; results returned by \code{\link{estimate_config}} run
#' with permuted datasets
#' @param thresholds vector of posterior probability thresholds for which to calculate eFDR
#' @param config numerical value; the column number of the configuration for which to
#' estimate eFDR values
#'
#' @return vector of the eFDR for each posterior probability threshold tested.
#'
#' @export
#'
efdr <- function(true_res,perm_res,thresholds=seq(0.99999,0.99900,-.00001),config=NULL){
  # last configuration (alternative for all data sources/types) is default for posterior prob. eFDR
  if (is.null(config)) config <- ncol(true_res$post_prob)
  pp <- true_res$post_prob[,config]
  pp0 <- sapply(perm_res,function(x) x$post_prob[,config])
  # calculate eFDR at each threshold
  efdr <- sapply(thresholds, function(x) mean(pp0>=x)/mean(pp>=x))
}


#' Posterior Probability threshold to control FDR
#'
#' Determine posterior probability threshold to use for controlling false discovery rate (FDR).
#' The function a vector of empirical FDR calculated for candidate posterior thresholds.
#' and returns a posterior probability threshold expected to control FDR at specified rate.
#'
#' @param efdr_res vector of eFDR values by potential posterior probability thresholds
#' @param fdr numerical value; desired/specified false discovery rate
#' @param tol the tolerance level for FDR control (i.e. half of the convergence window)
#' @param descending logical; denotes whether posterior probability thresholds are in
#' descending order
#' @param max_pp logical; denotes whether to return the maximum or minimum posterior probability
#' threshold which meets the FDR criteria
#'
#'
#' @return numerical value; posterior probability threshold if selected efdr was named
#' after probability value.
#' Otherwise, returns the index of the matching eFDR (with warning message).
#'
#' @export
#'
select_ppthresh <- function(efdr_res,fdr=0.05,tol=1e-3,descending=T,max_pp=T){
  all_matches <- which(efdr_res > (fdr-tol) & efdr_res < (fdr+tol))

  if (is.null(all_matches)) stop("No eFDR fuzzy-matched the specified FDR.\n
                             Expand posterior probability grid or relax tolerance threshold.")

  # if descending and max, or ascending and min, select maximum index
  if (descending==max_pp) select.idx <- max(all_matches)
  # otherwise, select minimum index
  else select.idx <- min(all_matches)

  # return index of matching threshold to control eFDR unless efdr_res named by threshold
  if (is.null(names(efdr_res)) | is.na(as.numeric(names(efdr_res)[select.idx]))){
    cat("Returning index of matching eFDR.")
    return(select.idx)
  } else return(as.numeric(names(efdr_res)[select.idx]))
}
