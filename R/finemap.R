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



#' Setup fine-mapping for a specified variant by index.
#'
#' For a specified variant (passed by index), set-up and run fine-mapping.
#' The function uses a data.table of leadSNPs to identify possible SNPs
#' to be conditioned on by the fine-mapping function, and determines which SNPs
#' will be conditioned on based on specified criteria. Returns the association
#' pattern with the highest posterior probability after fine-mapping.
#'
#' @param Primo_obj A list returned by running the \eqn{t}-statistic version
#' of Primo (i.e. \code{\link{Primo_tstat}})
#' @param IDs A data.table of the SNP and phenotype IDs corresponding to each row
#' of the Primo results stored in \code{Primo_obj}.
#' @param idx An integer of the index of the genetic variant (e.g. row of \code{Tstat_mod})
#' on which one wants to perform fine-mapping
#' @param leadSNPs_byRegion A data.table that stores the lead SNP of each phenotype
#' in each region of the Primo results. Also includes p-values for the lead SNPs.
#' See Details for format.
#' @param SNP_col A string of the column name of SNPs (must be "SNP" in current version).
#' @param pheno_cols A character vector of the column names of the phenotype ID columns.
#' @param snp.info A data.table reporting the chromosome and position of each SNP.
#' Columns should be: \code{SNP, CHR, POS}.
#' @param LD_mat A matrix of LD coefficients (\eqn{r^{2}}{r^2}).
#' Rows and columns should match the order of (\code{idx.snp}, \code{idx.leadsnps}).
#' @param LD_thresh A scalar corresponding to the LD coefficient (\eqn{r^{2}}{r^2})
#' threshold to be used for conditional analysis. Lead SNPs with \eqn{r^{2} <}{r^2 <}
#' \code{LD_thresh} with the \code{idx} variant will be conditioned on,
#' pending other criteria. Default value (1) signifies no consideration of LD in conditional analyses.
#' @param dist_thresh A scalar of the minimum number of base pairs away from the \code{idx} SNP
#' that a lead SNP must be in order to be considered for conditional analysis. Default value (0)
#' signifies no consideration of chromosomal distance in conditional analyses.
#' @param pval_thresh A scalar of the p-value threshold a lead SNP must be below
#' with the phenotype for which it is lead SNP in order to be considered for conditional analysis.
#' Default value (1) signifies no consideration of strength of effect in conditional analyses.
#' @param suffices A character vector of the suffices corresponding to columns in
#' \code{leadSNPs_byRegion}. See Details.
#'
#' @return The numerical value corresponding to the association pattern
#' with the highest posterior probability following fine-mapping adjustment.
#' The value returned, \eqn{k}, corresponds to the \eqn{k}-th column of \code{pis}
#' from the Primo output, and the \eqn{k}-th row of the \eqn{Q} matrix produced
#' by \code{\link{make_qmat}}.
#'
#' @details The following are additional details describing the input data.table
#' \code{leadSNPs_byRegion}. For \eqn{J} phenotypes being analyzed, the first
#' \eqn{J} columns of \code{leadSNPs_byRegion} should specify phenotype IDs.
#' Examples include gene symbol, CpG site name, trait name for GWAS, etc. The following
#' \eqn{2J} columns should hold the name and p-value of the lead SNP for each of the
#' \eqn{J} phenotypes. The column names should be of the form \code{leadSNP_x} and
#' \code{p-value_x}, where \code{x} is a suffix corresponding to the phenotype.
#'
#' @export
#'
fine_map_once <- function(Primo_obj,IDs,idx,leadSNPs_byRegion,SNP_col="SNP",pheno_cols,snp.info,LD_mat,LD_thresh=1,dist_thresh=0,pval_thresh=1,suffices=1:length(pheno_cols)){

  # IDs <- data.frame(IDs)
  # leadSNPs_byRegion <- data.frame(leadSNPs_byRegion)

  curr.IDs <- IDs[idx,]
  # curr.SNP <- curr.IDs[,get(SNP_col)]
  curr.SNP <- subset(curr.IDs,select=SNP_col)[[1]]
  # curr.SNP <- curr.IDs[,SNP_col]
  curr.Region <- merge(leadSNPs_byRegion,curr.IDs,by=pheno_cols)

  ## subset Primo results to the current region
  IDs_copy <- IDs
  IDs_copy$ObsNum <- 1:nrow(IDs_copy)
  IDs_copy <- merge(IDs_copy,curr.IDs,by=pheno_cols)
  curr_Region.idx <- IDs_copy$ObsNum
  Primo_obj <- subset_Primo_obj(Primo_obj,curr_Region.idx)
  IDs <- IDs[curr_Region.idx,]

  ## melt data so each lead SNP is its own row
  curr.Region_long <- melt(curr.Region,id.vars=pheno_cols,measure.vars = paste0("leadSNP_",suffices),
                                     value.name=SNP_col)
  curr.Region_long <- cbind(curr.Region_long,pval=melt(curr.Region,id.vars=pheno_cols,measure.vars = paste0("p-value_",suffices),
                                                                 value.name="pval")$pval)
  # curr.Region_long <- reshape2::melt(curr.Region,id.vars=pheno_cols,measure.vars = paste0("leadSNP_",suffices),
  #                          value.name=SNP_col)
  # curr.Region_long <- cbind(curr.Region_long,pval=reshape2::melt(curr.Region,id.vars=pheno_cols,measure.vars = paste0("p.value_",suffices),
  #                                                      value.name="pval")$pval)

  ## merge in chr and position to calculate distance between lead SNPs and SNP of interest
  # setkeyv(curr.Region_long,SNP_col)
  # setkeyv(snp.info,SNP_col)
  curr.Region_long <- merge(curr.Region_long,snp.info,by=SNP_col)
  curr.Region_long$dist <- abs(snp.info$POS[which(snp.info$SNP==curr.SNP)] - curr.Region_long$POS)
  ## merge in LD coefficients
  curr.Region_long$LD_r2 <- LD_mat[curr.SNP,subset(curr.Region_long,select=SNP_col)[[1]]]
  # curr.Region_long$LD_r2 <- LD_mat[curr.SNP,curr.Region_long[,SNP_col]]

  leadSNPs <- unique(subset(curr.Region_long, dist > dist_thresh & pval < pval_thresh & LD_r2 < LD_thresh)$SNP)
  # leadSNPs <- unique(subset(curr.Region_long, dist > dist_thresh & pval < pval_thresh & LD_r2 < LD_thresh)[,SNP_col])

  ## index of SNP of interest
  # idx.snp <- which(IDs[,get(SNP_col)]==curr.SNP)
  idx.snp <- which(subset(IDs,select=SNP_col)[[1]]==curr.SNP)
  # idx.snp <- which(IDs[,SNP_col]==curr.SNP)

  if(length(leadSNPs)==0){
    return(which.max(Primo_obj$post_prob[idx.snp,]))
  } else{

    ## get snp_indices for other snps to adjust for
    idx.leadsnps <- NULL
    for(j in 1:length(leadSNPs)){
      # idx.leadsnps <- c(idx.leadsnps, which(IDs[,get(SNP_col)]==leadSNPs[j]))
      idx.leadsnps <- c(idx.leadsnps, which(subset(IDs,select=SNP_col)[[1]]==leadSNPs[j]))
      # idx.leadsnps <- c(idx.leadsnps, which(IDs[,SNP_col]==leadSNPs[j]))
    }

    ## run fine-mapping
    sp <- fine_map(idx.snp,idx.leadsnps,LD_mat[c(curr.SNP,leadSNPs),c(curr.SNP,leadSNPs)],Primo_obj)
    return(sp)
  }
}
