#' @aliases Primo-package
#'
#' @details
#' To estimate the posterior probabilities of belonging to each association pattern,
#' either user the wrapper function \code{\link{Primo}}, or directly call
#' \code{\link{Primo_tstat}} (for \eqn{t}-statistics) or \code{\link{Primo_pval}} (for \eqn{P}-values).
#'
#'
"_PACKAGE"

#' @useDynLib Primo, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom magrittr %>%
#' @exportPattern "^[[:alpha:]]+"
NULL
