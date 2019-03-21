#' Perform conditional analysis for a specified variant.
#'
#' For a specified variant, re-estimate the posterior probabilities of association
#' patterns, conditioning on other specified variants (e.g. lead SNPs for phenotypes
#' being studied). Returns the association pattern with the highest
#' posterior probability after conditional analysis.
#'
#' @param idx_snp integer matching the index of the genetic variant (e.g. row of \code{Tstat_mod})
#' for which to perform conditional analysis.
#' @param idx_leadsnps vector of indices of the leading snps (e.g. rows of \code{Tstat_mod})
#' on which to condition.
#' @param LD_mat matrix of LD coefficients (\eqn{r^{2}}{r^2}).
#' Rows and columns should match the order of (\code{idx_snp}, \code{idx_leadsnps}).
#' @param Primo_obj list returned by running the \eqn{t}-statistic version
#' of Primo (i.e. \code{\link{Primo_tstat}} or \code{\link{Primo_modT}})
#'
#' @return The integer corresponding to the association pattern
#' with the highest posterior probability following conditional analysis.
#' The value returned, \eqn{k}, corresponds to the \eqn{k}-th column of \code{pis}
#' from the Primo output, and the \eqn{k}-th row of the \eqn{Q} matrix produced
#' by \code{\link{make_qmat}}.
#'
#'
#' @export
#'
Primo_conditional <- function(idx_snp,idx_leadsnps,LD_mat,Primo_obj){
  n_leadsnps<-length(idx_leadsnps)
  zi<-Primo_obj$Tstat_mod[c(idx_snp,idx_leadsnps),]
  d<-ncol(zi)
  Q <- make_qmat(1:d)
  V <- Primo_obj$V_mat[c(idx_snp,idx_leadsnps),]
  mdf_sd<-Primo_obj$mdf_sd_mat[c(idx_snp,idx_leadsnps),]
  post_prob<-Primo_obj$post_prob
  pis<-Primo_obj$pis
  v2<-NULL
  Gamma<-Primo_obj$Gamma
  for(i in 1:n_leadsnps){
    q2<-Q[which.max(post_prob[idx_leadsnps[i],]),]
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

  return(sp)
}



#' Set up conditional analysis for a specified variant by index.
#'
#' For a specified variant (passed by index), set-up and run conditional analysis.
#' The function uses a data.frame of lead SNPs to identify possible SNPs
#' for conditional analysis, and determines which SNPs
#' will be conditioned on based on specified criteria. Returns the association
#' pattern with the highest posterior probability after conditional analysis.
#'
#' @param Primo_obj list returned by running the \eqn{t}-statistic version
#' of Primo (i.e. \code{\link{Primo_tstat}})
#' @param IDs data.frame of the SNP and phenotype IDs corresponding to each row
#' of \code{Primo_obj}.
#' @param idx integer specifying the index of the genetic variant (e.g. row of \code{Tstat_mod})
#' on which one wants to perform conditional analysis.
#' @param leadsnps_region data.frame detailing the lead SNP of each phenotype
#' in each region of \code{Primo_obj} and associated \eqn{P}-values for the lead SNPs.
#' See Details.
#' @param snp_col string of the column name of SNPs/variants.
#' @param pheno_cols character vector of the column names of the phenotype ID columns.
#' @param snp_info data.frame reporting the chromosome and position of each SNP.
#' Columns should be: \code{SNP, CHR, POS}.
#' @param LD_mat matrix of LD coefficients (\eqn{r^{2}}{r^2}). Row and column names
#' should be SNP/variant names (e.g. in \code{snp_col}).
#' @param LD_thresh scalar corresponding to the LD coefficient (\eqn{r^{2}}{r^2})
#' threshold to be used for conditional analysis. Lead SNPs with \eqn{r^{2} <}{r^2 <}
#' \code{LD_thresh} with the \code{idx} variant will be conditioned on.
#' Default value (1) signifies no consideration of LD in conditional analyses.
#' @param dist_thresh scalar of the minimum number of base pairs away from the \code{idx} SNP
#' that a lead SNP must be in order to be conditioned on. Default value (0)
#' signifies no consideration of chromosomal distance in conditional analyses.
#' @param pval_thresh scalar of the \eqn{P}-value threshold a lead SNP must be below
#' with the phenotype for which it is lead SNP in order to be conditionaed on.
#' Default value (1) signifies no consideration of strength of effect in conditional analyses.
#' @param suffices character vector of the suffices corresponding to columns in
#' \code{leadsnps_region}. See Details.
#'
#' @inherit Primo_conditional return
#'
#' @details The following are additional details describing the input
#' \code{leadsnps_region}. For \eqn{J} phenotypes being analyzed, the first
#' \eqn{J} columns of \code{leadsnps_region} should specify phenotype IDs.
#' Examples include gene symbol, CpG site name, trait name for GWAS, etc. The following
#' \eqn{2*J} columns should hold the name and p-value of the lead SNP for each of the
#' \eqn{J} phenotypes. The column names should be of the form \code{leadSNP_x} and
#' \code{pvalue_x}, where \code{x} is a suffix corresponding to the phenotype.
#'
#' @export
#'
run_conditional <- function(Primo_obj,IDs,idx,leadsnps_region,snp_col="SNP",pheno_cols,snp_info,LD_mat,LD_thresh=1,dist_thresh=0,pval_thresh=1,suffices=1:length(pheno_cols)){

  curr.IDs <- IDs[idx,]
  curr.SNP <- curr.IDs[,snp_col]
  curr.Region <- merge(leadsnps_region,curr.IDs,by=pheno_cols)

  ## subset Primo results to the current region
  IDs_copy <- IDs
  IDs_copy$ObsNum <- 1:nrow(IDs_copy)
  IDs_copy <- merge(IDs_copy,curr.IDs,by=pheno_cols)
  curr_Region.idx <- IDs_copy$ObsNum
  Primo_obj <- subset_Primo_obj(Primo_obj,curr_Region.idx)
  IDs <- IDs[curr_Region.idx,]

  ## melt data so each lead SNP is its own row
  curr.Region_long <- reshape2::melt(curr.Region,id.vars=pheno_cols,measure.vars = paste0("leadSNP_",suffices),
                           value.name=snp_col)
  curr.Region_long <- cbind(curr.Region_long,pval=reshape2::melt(curr.Region,id.vars=pheno_cols,measure.vars = paste0("pvalue_",suffices),
                                                       value.name="pval")$pval)

  ## merge in chr and position to calculate distance between lead SNPs and SNP of interest
  curr.Region_long <- merge(curr.Region_long,snp_info,by=snp_col)
  curr.Region_long$dist <- abs(snp_info$POS[which(snp_info$SNP==curr.SNP)] - curr.Region_long$POS)
  ## merge in LD coefficients
  curr.Region_long$LD_r2 <- LD_mat[curr.SNP,curr.Region_long[,snp_col]]

  # leadSNPs <- unique(subset(curr.Region_long, dist > dist_thresh & pval <= pval_thresh & LD_r2 < LD_thresh)[,snp_col])
  keep_idx <- which(curr.Region_long$dist > dist_thresh & curr.Region_long$pval <= pval_thresh & curr.Region_long$LD_r2 < LD_thresh)
  leadSNPs <- unique(curr.Region_long[keep_idx,snp_col])

  ## index of SNP of interest
  idx_snp <- which(IDs[,snp_col]==curr.SNP)

  if(length(leadSNPs)==0){
    return(which.max(Primo_obj$post_prob[idx_snp,]))
  } else{

    ## get snp_indices for other snps to adjust for
    idx_leadsnps <- NULL
    for(j in 1:length(leadSNPs)){
      idx_leadsnps <- c(idx_leadsnps, which(IDs[,snp_col]==leadSNPs[j]))
    }

    ## run fine-mapping
    sp <- Primo::Primo_conditional(idx_snp,idx_leadsnps,LD_mat[c(curr.SNP,leadSNPs),c(curr.SNP,leadSNPs)],Primo_obj)
    return(sp)
  }
}

#' Setup fine-mapping for a specified variant by index, using data.tables.
#'
#' For a specified variant (passed by index), set-up and run conditional analysis.
#' The function uses a data.table of lead SNPs to identify possible SNPs
#' for conditional analysis, and determines which SNPs
#' will be conditioned on based on specified criteria. Returns the association
#' pattern with the highest posterior probability after conditional analysis.
#' The version utilizes data.tables for \code{IDs}, \code{leadsnps_region}
#' and \code{snp_info} for faster processing.
#'
#' @param Primo_obj list returned by running the \eqn{t}-statistic version
#' of Primo (i.e. \code{\link{Primo_tstat}})
#' @param IDs data.table of the SNP and phenotype IDs corresponding to each row
#' of the Primo results stored in \code{Primo_obj}.
#' @param idx integer of the index of the genetic variant (e.g. row of \code{Tstat_mod})
#' on which one wants to perform fine-mapping
#' @param leadsnps_region data.table that stores the lead SNP of each phenotype
#' in each region of the Primo results. Also includes p-values for the lead SNPs.
#' See Details for format.
#' @param snp_col string of the column name of SNPs/variants.
#' @param pheno_cols character vector of the column names of the phenotype ID columns.
#' @param snp_info data.table reporting the chromosome and position of each SNP.
#' Columns should be: \code{SNP, CHR, POS}.
#' @param LD_mat matrix of LD coefficients (\eqn{r^{2}}{r^2}). Row and column names
#' should be SNP/variant names (e.g. in \code{snp_col}).
#' @param LD_thresh scalar corresponding to the LD coefficient (\eqn{r^{2}}{r^2})
#' threshold to be used for conditional analysis. Lead SNPs with \eqn{r^{2} <}{r^2 <}
#' \code{LD_thresh} with the \code{idx} variant will be conditioned on.
#' Default value (1) signifies no consideration of LD in conditional analyses.
#' @param dist_thresh scalar of the minimum number of base pairs away from the \code{idx} SNP
#' that a lead SNP must be in order to be conditioned on. Default value (0)
#' signifies no consideration of chromosomal distance in conditional analyses.
#' @param pval_thresh scalar of the \eqn{P}-value threshold a lead SNP must be below
#' with the phenotype for which it is lead SNP in order to be conditionaed on.
#' Default value (1) signifies no consideration of strength of effect in conditional analyses.
#' @param suffices character vector of the suffices corresponding to columns in
#' \code{leadsnps_region}. See Details.
#'
#' @inherit Primo_conditional return
#'
#' @inherit run_conditional details
#'
#' @export
#'
run_conditional_dt <- function(Primo_obj,IDs,idx,leadsnps_region,snp_col="SNP",pheno_cols,snp_info,LD_mat,LD_thresh=1,dist_thresh=0,pval_thresh=1,suffices=1:length(pheno_cols)){

  base::requireNamespace("data.table")

  curr.IDs <- IDs[idx,]
  # curr.SNP <- curr.IDs[,get(snp_col)]
  curr.SNP <- subset(curr.IDs,select=snp_col)[[1]]
  curr.Region <- merge(leadsnps_region,curr.IDs,by=pheno_cols)

  ## subset Primo results to the current region
  IDs_copy <- IDs
  IDs_copy$ObsNum <- 1:nrow(IDs_copy)
  IDs_copy <- merge(IDs_copy,curr.IDs,by=pheno_cols)
  curr_Region.idx <- IDs_copy$ObsNum
  Primo_obj <- subset_Primo_obj(Primo_obj,curr_Region.idx)
  IDs <- IDs[curr_Region.idx,]

  ## melt data so each lead SNP is its own row
  curr.Region_long <- data.table::melt(curr.Region,id.vars=pheno_cols,measure.vars = paste0("leadSNP_",suffices),
                                     value.name=snp_col)
  curr.Region_long <- cbind(curr.Region_long,pval=data.table::melt(curr.Region,id.vars=pheno_cols,measure.vars = paste0("pvalue_",suffices),
                                                                 value.name="pval")$pval)

  ## merge in chr and position to calculate distance between lead SNPs and SNP of interest
  data.table::setkeyv(curr.Region_long,snp_col)
  data.table::setkeyv(snp_info,snp_col)
  curr.Region_long <- merge(curr.Region_long,snp_info)
  curr.Region_long$dist <- abs(snp_info$POS[which(snp_info$SNP==curr.SNP)] - curr.Region_long$POS)
  ## merge in LD coefficients
  curr.Region_long$LD_r2 <- LD_mat[curr.SNP,subset(curr.Region_long,select=snp_col)[[1]]]

  # leadSNPs <- unique(subset(curr.Region_long, dist > dist_thresh & pval <= pval_thresh & LD_r2 < LD_thresh,select=snp_col)[[1]])
  keep_idx <- which(curr.Region_long$dist > dist_thresh & curr.Region_long$pval <= pval_thresh & curr.Region_long$LD_r2 < LD_thresh)
  # leadSNPs <- unique(curr.Region_long[keep_idx,..snp_col][[1]])
  # leadSNPs <- unique(subset(curr.Region_long,select=snp_col)[[1]])
  leadSNPs <- unique(curr.Region_long$SNP[keep_idx])

  ## index of SNP of interest
  # idx_snp <- which(IDs[,get(snp_col)]==curr.SNP)
  idx_snp <- which(subset(IDs,select=snp_col)[[1]]==curr.SNP)
  # idx_snp <- which(IDs[,..snp_col][[1]]==curr.SNP)

  if(length(leadSNPs)==0){
    return(which.max(Primo_obj$post_prob[idx_snp,]))
  } else{

    ## get snp_indices for other snps to adjust for
    idx_leadsnps <- NULL
    for(j in 1:length(leadSNPs)){
      # idx_leadsnps <- c(idx_leadsnps, which(IDs[,get(snp_col)]==leadSNPs[j]))
      idx_leadsnps <- c(idx_leadsnps, which(subset(IDs,select=snp_col)[[1]]==leadSNPs[j]))
      # idx_leadsnps <- c(idx_leadsnps, which(IDs[,..snp_col][[1]]==leadSNPs[j]))
    }

    ## run fine-mapping
    sp <- Primo::Primo_conditional(idx_snp,idx_leadsnps,LD_mat[c(curr.SNP,leadSNPs),c(curr.SNP,leadSNPs)],Primo_obj)
    return(sp)
  }
}

