#' @title PRS_tf
#' @description  This function is used to generate the GWAS data for TL-Multi which were applied to conduct multiethnic PRS prediction
#' @details This function will return the multiethnic PRS calculated by TL-Multi
#' @param v_e The validation result of Lassosum for the auxiliary informative population
#' @param out_e The Lassosum result of the auxiliary informative population
#' @param ss_a The GWAS data for the target population
#' @param XTX The LD region of the target population
#' @return A list of multiethnic PRS calculated by TL-Multi.
#'
#'
#'
#'snp = snp,

ss_tl <- function(v_e, out_e, ss_a, XTX = XTX){

  
  beta_eur = v_e$best.beta # best beta from validation
  beta_eur = data.frame(beta_eur, out_e$sumstats$pos)
  
  # Obtain Asian summary statistics from input
  df_tl = ss_a
  names(df_tl)[names(df_tl) == "beta"] = "beta_asn"
  
  beta_eur = data.frame(beta_eur, out_e$sumstats$pos, 
                        out_e$sumstats$chr)
  
  
  # Merge the dataframes by their positions
  df_tl = merge(df_tl, beta_eur, by.x=c("bp",'chr'), 
                by.y=c("out_e.sumstats.pos",'out_e.sumstats.chr'),  sort=F)
  
  best_beta <- data.frame(v_e$best.beta, out_e$sumstats$pos, 
                          out_e$sumstats$chr)
  colnames(best_beta) <- c('beta_e', 'bp', 'chr')
  
  ss_tl <-  merge(df_tl, snplist, by.x="rsid", by.y="snplist", sort=F)
  ZTZ <-  XTX[ss_tl$rsid, ss_tl$rsid]

  ss_tl$beta_tl <-  ss_tl$beta_a - as.vector(ss_tl$beta_e %*% ZTZ)
  ss_tl$z <-  ss_tl$beta_tl / ss_tl$se
  ss_tl$pvalue <-  2 * pnorm(-abs(ss_tl$z))





  return(ss_tl)

}


