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

#' Set up conditional analysis for specified variant(s) by index.
#'
#' For specified variant(s) (passed by index), set-up and run conditional analysis.
#' The function uses a data.table of lead SNPs to identify possible SNPs
#' for conditional analysis, and determines which SNPs
#' will be conditioned on based on specified criteria. Returns the association
#' pattern with the highest posterior probability for each specified variant
#' after conditional analysis.
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
#' Columns must include: \code{SNP, CHR, POS}.
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
#' with the phenotype for which it is lead SNP in order to be conditioned on.
#' Default value (1) signifies no consideration of strength of effect in conditional analyses.
#' @param suffices character vector of the suffices corresponding to columns in
#' \code{leadsnps_region}. See Details.
#'
#' @inherit Primo_conditional return
#'
#' @return An integer vector corresponding to the association pattern
#' with the highest posterior probabilities for each SNP variant
#' represented by \code{idx}, following conditional analysis.
#' Each value returned, \eqn{k}, corresponds to the \eqn{k}-th column of \code{pis}
#' from the Primo output, and the \eqn{k}-th row of the \eqn{Q} matrix produced
#' by \code{\link{make_qmat}}.
#'
#' @export
#'
run_conditional <- function(Primo_obj,IDs,idx,leadsnps_region,snp_col="SNP",pheno_cols,snp_info,LD_mat,LD_thresh=1,dist_thresh=0,pval_thresh=1,suffices=1:length(pheno_cols)){

  base::requireNamespace("data.table")

  if(!data.table::is.data.table(IDs)) IDs <- data.table::data.table(IDs)
  if(!data.table::is.data.table(leadsnps_region)) IDs <- data.table::data.table(leadsnps_region)
  if(!data.table::is.data.table(snp_info)) IDs <- data.table::data.table(snp_info)

  sp_vec <- NULL

  for(i in idx){
    ## get current region
    curr.IDs <- IDs[i,]
    curr.SNP <- subset(curr.IDs,select=snp_col)[[1]]
    curr.Region <- merge(leadsnps_region,curr.IDs,by=pheno_cols)

    ## subset Primo results to the current region
    IDs_copy <- IDs
    IDs_copy$ObsNum <- 1:nrow(IDs_copy)
    IDs_copy <- merge(IDs_copy,curr.IDs,by=pheno_cols)
    curr_Region.idx <- IDs_copy$ObsNum
    Primo_obj_sub <- subset_Primo_obj(Primo_obj,curr_Region.idx)
    IDs_copy <- IDs[curr_Region.idx,]

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

    ## determine which lead SNPs meet criteria
    keep_idx <- which(curr.Region_long$dist > dist_thresh & curr.Region_long$pval <= pval_thresh & curr.Region_long$LD_r2 < LD_thresh)
    leadSNPs <- unique(subset(curr.Region_long[keep_idx,],select=snp_col)[[1]])

    ## index of SNP of interest
    idx_snp <- which(subset(IDs_copy,select=snp_col)[[1]]==curr.SNP)

    if(length(leadSNPs)==0){
      if(is.matrix(Primo_obj_sub$post_prob)){
        sp_vec <- c(sp_vec,which.max(Primo_obj_sub$post_prob[idx_snp,]))
      } else sp_vec <- c(sp_vec,which.max(Primo_obj_sub$post_prob))
    } else{

      ## get snp_indices for other snps to adjust for
      idx_leadsnps <- NULL
      for(j in 1:length(leadSNPs)){
        idx_leadsnps <- c(idx_leadsnps, which(subset(IDs_copy,select=snp_col)[[1]]==leadSNPs[j]))
      }

      ## run conditional analysis
      sp <- Primo::Primo_conditional(idx_snp,idx_leadsnps,LD_mat[c(curr.SNP,leadSNPs),c(curr.SNP,leadSNPs)],Primo_obj_sub)
      sp_vec <- c(sp_vec,sp)
    }
  }

  return(sp_vec)
}

#' Set up conditional analysis for known complex trait-associated variants.
#'
#' For specified, known complex trait-associated (GWAS) variant(s),
#' set-up and run conditional analysis.
#' The function identifies lead omics SNPs to consider for conditional analysis,
#' and determines which SNPs will be conditioned on for each GWAS variant
#' based on specified criteria. Returns a data.frame with posterior probabilities
#' for collapsed association patterns and results from conditional analysis,
#' as well as estimated FDR for each collapsed association pattern at a
#' specified posterior probability threshold.
#'
#' @param Primo_obj list returned by running the \eqn{t}-statistic version
#' of Primo (i.e. \code{\link{Primo_tstat}})
#' @param IDs data.frame of the SNP and phenotype IDs corresponding to each row
#' of the Primo results stored in \code{Primo_obj}.
#' @param gwas_snps character vector of known trait-associated (GWAS) SNPs.
#' @param pvals matrix of \eqn{P}-values from test statistics.
#' @param LD_mat matrix of LD coefficients (\eqn{r^{2}}{r^2}). Row and column names
#' should be SNP/variant names (i.e matching those present in \code{IDs}).
#' @param snp_info data.frame reporting the chromosome and position of each SNP.
#' Columns must include: \code{SNP, CHR, POS}.
#' @param pp_thresh scalar of the posterior probability threshold used for significance.
#' @param LD_thresh scalar corresponding to the LD coefficient (\eqn{r^{2}}{r^2})
#' threshold to be used for conditional analysis. Lead omics SNPs with \eqn{r^{2} <}{r^2 <}
#' \code{LD_thresh} with the GWAS SNP will be conditioned on.
#' @param dist_thresh scalar of the minimum number of base pairs away from the GWAS SNP
#' that a lead SNP must be in order to be conditioned on.
#' @param pval_thresh scalar of the \eqn{P}-value threshold a lead SNP must be below
#' with the phenotype for which it is lead SNP in order to be conditioned on.
#'
#'
#' @return A list with two elements, \code{pp_grouped} and \code{fdr}.
#'
#' \code{fdr} is a named vector of the estimated false discovery rates (FDR)
#' for each collapsed association pattern at the posterior probability
#' threshold, \code{pp_thresh}.
#'
#' \code{pp_grouped} is a data.frame with the following information:
#'
#' \itemize{
#'   \item SNP and trait identifiers corresponding to each observation
#'   \item posterior probabilities of the collapsed association patterns
#'   ("GWAS + at least x omics trait(s)")
#'   \item number of omics traits with which the SNP was associated before conditional analysis
#'   (at posterior probability \code{> pp_thresh})
#'   \item number of omics traits the SNP is associated with after conditional analysis
#'   \item the top association pattern after conditional analysis
#' }
#'
#' @export
#'
run_conditional_gwas <- function(Primo_obj,IDs,gwas_snps,pvals,LD_mat,snp_info,pp_thresh,LD_thresh=0.9,dist_thresh=5e3,pval_thresh=1e-3){

  snp_col <- colnames(IDs)[1]
  pheno_cols <- colnames(IDs)[2:ncol(IDs)]

  # gwas_idx <- which(IDs[,snp_col] %in% gwas_snps)
  ## use subset command to allow IDs to be either data.table or data.frame
  gwas_idx <- which(subset(IDs, select=snp_col)[[1]] %in% gwas_snps)

  colnames(pvals) <- paste0("pvalue_",1:ncol(pvals))

  ## determine lead omics SNPs
  ID_pvals <- data.table::data.table(IDs,pvals)
  leadsnps_region <- Primo::find_leadsnps(data=ID_pvals,snp_col=snp_col,pheno_cols=pheno_cols,
                                          stat_cols=colnames(ID_pvals)[(ncol(IDs)+1):ncol(ID_pvals)],data_type="pvalue")

  ## determine top pattern after conditional analysis
  sp_vec <- Primo::run_conditional(Primo_obj,IDs,idx=gwas_idx,leadsnps_region,snp_col,
                                   pheno_cols,snp_info,LD_mat,LD_thresh,dist_thresh,pval_thresh)

  ## collapsed probabilities
  PP_grouped <- Primo::collapse_pp_num(Primo_obj$post_prob[gwas_idx,],req_idx=1,prefix="pp_nQTL_ge")
  PP_grouped <- PP_grouped[,-1] ## drop 0qtl+GWAS since that is not goal of analysis

  ## number of QTLs, per collapsed probabilities
  PP_grouped <- cbind(PP_grouped,nQTL_orig=0)
  for(k in 1:(ncol(pvals)-1)){
    PP_grouped[,"nQTL_orig"] <- PP_grouped[,"nQTL_orig"] + as.numeric(PP_grouped[,paste0("pp_nQTL_ge",k)] > pp_thresh)
  }

  ## determine the number of QTLs in each corresponding pattern
  myQ <- Primo::make_qmat(1:ncol(pvals))
  nQTL_byQrow <- rowSums(myQ[,2:ncol(myQ)])
  nQTL_final <- nQTL_byQrow[sp_vec]

  ## add in sp_vec and final nQTL
  PP_grouped <- cbind(PP_grouped,nQTL_final=nQTL_final,top_pattern=sp_vec)
  PP_grouped[,"nQTL_final"] <- pmin(PP_grouped[,"nQTL_orig"],PP_grouped[,"nQTL_final"])

  ## calculate estimated FDR
  fdr <- NULL
  for(i in 1:(ncol(pvals)-1)){
    pp_col <- paste0("pp_nQTL_ge",i)
    next_fdr <- Primo::calc_fdr_conditional(PP_grouped[,pp_col],thresh=pp_thresh,
                                            fail_idx=which(PP_grouped[,"nQTL_orig"]>=i & PP_grouped[,"nQTL_final"] < i))
    fdr <- c(fdr,next_fdr)
  }
  names(fdr) <- paste0("nQTL_ge",1:(ncol(pvals)-1))

  PP_grouped <- data.frame(IDs[gwas_idx,],PP_grouped,stringsAsFactors = F)

  PP_grouped <- cbind(IDs[gwas_idx,], PP_grouped)

  return(list(pp_grouped=PP_grouped,fdr=fdr))
}
