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


#' Find the lead SNP for each phenotype in each region.
#'
#' Subset results from Primo output based on a vector of indices.
#'
#' @param data A data.table. Each row will be a SNP-phenotype combination
#' with statistics necessary to determine the lead SNP in each phenotype region.
#' @param SNP_col Character string of the column name of the SNP (must be "SNP"
#' for current version).
#' @param pheno_cols Character vector of the column names of the phenotypes.
#' @param stat_cols Character vector of the column names of statistics to be
#' used to determine lead SNPs.
#' @param data_type Character string denoting type of statistics being used. Must be
#' either "pvalue" or "tstat".
#' @param suffices A character vector denoting suffices to use for the names of the
#' lead SNP columns (optional). If \code{NULL}, consecutive integers will be assigned.
#'
#' @return A data.table containing information about the lead SNPs and
#' associated statistics. The columns will be \code{pheno_cols} followed by two columns
#' for each phenotype: the name of the lead SNP and the value of the statistic for
#' that lead SNP in that phenotype.
#'
#' @export
#'
find_leadSNPs <- function(data,SNP_col,pheno_cols,stat_cols,data_type="pvalue",suffices=NULL){

  require(data.table)

  cat("class: ", class(data))

  setkeyv(data,pheno_cols)

  if(is.null(suffices)) suffices <- 1:length(stat_cols)

  if(data_type=="pvalue"){

    ## find SNP with minimum p-value for the first phenotype in the region
    leadSNPs_byRegion <- data[data[,.I[get(stat_cols[1])==min(get(stat_cols[1]))],by=key(data)]$V1]
    leadSNPs_byRegion<- subset(leadSNPs_byRegion, select=c(pheno_cols,"SNP",stat_cols[1]))
    colnames(leadSNPs_byRegion)[ncol(leadSNPs_byRegion)-1] <- paste0("leadSNP_",suffices[1])

    ## drop duplicates if they exist (cases of perfect LD)
    leadSNPs_byRegion <- leadSNPs_byRegion[!duplicated(leadSNPs_byRegion,by=key(leadSNPs_byRegion))]

    for(i in 2:length(stat_cols)){

      ## find SNP with minimum p-value for the current phenotype in the region
      topSNP_currPheno <- data[data[,.I[get(stat_cols[i])==min(get(stat_cols[i]))],by=key(data)]$V1]
      topSNP_currPheno<- subset(topSNP_currPheno, select=c(pheno_cols,"SNP",stat_cols[i]))
      colnames(topSNP_currPheno)[ncol(topSNP_currPheno)-1] <- paste0("leadSNP_",suffices[i])

      ## drop duplicates if they exist (cases of perfect LD) and merge
      topSNP_currPheno <- topSNP_currPheno[!duplicated(topSNP_currPheno,by=key(topSNP_currPheno))]
      leadSNPs_byRegion <- merge(leadSNPs_byRegion,topSNP_currPheno)
    }

  } else if(data_type=="tstat"){

    ## find SNP with maximum abs(t-statistic) for the first phenotype in the region
    leadSNPs_byRegion <- data[data[,.I[get(stat_cols[1])==max(abs(get(stat_cols[1])))],by=key(data)]$V1]
    leadSNPs_byRegion<- subset(leadSNPs_byRegion, select=c(pheno_cols,"SNP",stat_cols[1]))
    colnames(leadSNPs_byRegion)[ncol(leadSNPs_byRegion)-1] <- paste0("leadSNP_",suffices[1])

    ## drop duplicates if they exist (cases of perfect LD)
    leadSNPs_byRegion <- leadSNPs_byRegion[!duplicated(leadSNPs_byRegion,by=key(leadSNPs_byRegion))]

    for(i in 2:length(stat_cols)){

      ## find SNP with maximum abs(t-statistic) for the current phenotype in the region
      topSNP_currPheno <- data[data[,.I[get(stat_cols[i])==min(get(stat_cols[i]))],by=key(data)]$V1]
      topSNP_currPheno<- subset(topSNP_currPheno, select=c(pheno_cols,"SNP",stat_cols[i]))
      colnames(topSNP_currPheno)[ncol(topSNP_currPheno)-1] <- paste0("leadSNP_",suffices[i])

      ## drop duplicates if they exist (cases of perfect LD)
      topSNP_currPheno <- topSNP_currPheno[!duplicated(topSNP_currPheno,by=key(topSNP_currPheno))]
      leadSNPs_byRegion <- merge(leadSNPs_byRegion,topSNP_currPheno)
    }

  } else{

    stop("data_type must be either 'pvalue' or 'tstat'.")

  }


  return(leadSNPs_byRegion)
}
