#' @details
#' To estimate the posterior probabilities of belonging to each association pattern,
#' either user the wrapper function \code{\link{Primo}}, or directly call
#' \code{\link{Primo_tstat}} (for t-statistics) or \code{\link{Primo_pval}} (for p-values).
#'
"_PACKAGE"

#' @useDynLib primo, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL
