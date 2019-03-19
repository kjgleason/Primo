#' Calculate the empirical false discovery rate (FDR).
#'
#' Calculate the empirical false discovery rate (FDR) given
#' a vector of posterior probabilities and a specified threshold.
#'
#' @param post_prob numeric vector of posterior probabilities.
#' @param threshold numeric value of the posterior probability threshold at which to
#' calculate the empirical FDR.
#'
#' @return A numeric value of the empirical false discovery rate (FDR).
#'
#' @export
#'
calc_fdr <- function(post_prob,threshold){
  post_prob <- post_prob[which(post_prob > threshold)]
  return(mean(1-post_prob))
}


#' Calculate the empirical false discovery rate (FDR) after fine-mapping.
#'
#' Calculate the empirical false discovery rate (FDR) given
#' a vector of posterior probabilities and a specified threshold,
#' adjusting for false positives identified by conditional analysis.
#'
#' @param post_prob numeric vector of posterior probabilities.
#' @param threshold numerical value of the posterior probability threshold at which to
#' calculate the empirical FDR.
#' @param fail_idx integer vector of the indices of observations which "failed"
#' conditional analysis. These are the observations where the association pattern with the
#' highest posterior probability changed after conditioning. See important Note.
#'
#' @inherit calc_fdr return
#'
#' @section Note: Note that if using collapsed posterior probability categories,
#' an observation where the association pattern with the highest posterior probability
#' changes may not indicate "failure" if the new pattern still fits the
#' description of the collapsed category. For example, if testing "associated
#' with at least 1 trait," an observation that changes from "associated with 2"
#' to "associated with 1" may not be considered failure. In such cases, the user
#' may wish \code{fail_idx} to represent cases where the highest
#' association pattern no longer fits the description of the collapsed category).
#'
#' @export
#'
calc_fdr_conditional <- function(post_prob,threshold,fail_idx){

  ## store indices for observations that would have been fine-mapped according to threshold
  conditioned_idx <- which(post_prob > threshold)

  ## set PP=0 for those that failed fine-mapping (considered false positives)
  post_prob[fail_idx] <- 0

  post_prob <- post_prob[conditioned_idx]
  return(mean(1-post_prob))
}

#' Calculate the empirical false discovery rate (FDR) for a matrix column.
#'
#' Calculate the empirical false discovery rate (FDR) of a
#' column in a matrix of posterior probabilities given a specified threshold.
#'
#' @param post_prob numeric matrix of posterior probabilities.
#' @param threshold numeric value of the posterior probability threshold at which to
#' calculate the empirical FDR.
#' @param col_name character string of the column name.
#' @param col_idx integer of the column index.
#'
#' @inherit calc_fdr return
#'
#' @details At least one of \code{col_name} or \code{col_idx} must be specified.
#'
#' @export
#'
calc_fdr_col <- function(post_prob,threshold,col_name=NULL,col_idx=NULL){

  if(!is.null(col_name)){

    return(calc_fdr(post_prob[,"col_name"],threshold))

  } else if(!is.null(col_idx)){

    return(calc_fdr(post_prob[,col_idx],threshold))

  } else{

    stop("At least one of col_name or col_idx must be specified (non-NULL).")

  }
}

#' Calculate the empirical false discovery rate (FDR) for each matrix column.
#'
#' Calculate the empirical false discovery rate (FDR) of specified
#' columns in a matrix of posterior probabilities given a specified threshold.
#'
#' @param post_prob numeric matrix of posterior probabilities.
#' @param threshold numeric value of the posterior probability threshold at which to
#' calculate the empirical FDR.
#' @param col_names character vector of column names.
#' @param col_idxs integer vector of column indices.
#'
#' @inherit calc_fdr return
#'
#' @details At least one of \code{col_names} or \code{col_idxs} must be specified.
#'
#' @export
#'
calc_fdr_multi <- function(post_prob,threshold,col_names=NULL,col_idxs=NULL){

  if(is.null(col_names) & is.null(col_idxs)) stop("At least one of col_names or col_idxs must be specified (non-NULL).")

  if(!is.null(col_names)){

    fdr <- NULL
    for(i in 1:length(col_names)){
      fdr <- c(fdr,calc_fdr_oneCol(post_prob,threshold,col_name=col_names[i]))
    }

    names(fdr) <- col_names

  } else {

    fdr <- NULL
    for(i in 1:length(col_idxs)){
      fdr <- c(fdr,calc_fdr_oneCol(post_prob,threshold,col_idx=col_idxs[i]))
    }

  }

  return(fdr)
}
