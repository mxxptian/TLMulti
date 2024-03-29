
#' @title prepare_data
#' @description  This function is used to prepare data for generating the GWAS
#'   data for TL-Multi. It was based on the lassosum proposed by Mak et al
#'   (2017).
#' @param ss_e The GWAS data for the informative auxiliary population
#' @param ss_a The GWAS data for the target population
#' @param pheno The phenotype of the target population
#' @param ref.bfile.aux The reference panel for the informative auxiliary
#'   population
#' @param test.bfile.aux The testing panel for the informative auxiliary
#'   population
#' @param sample_t the test samples from data
#' @param ref.bfile.tar The reference panel for the target population
#' @param test.bfile.tar The testing panel for the target population
#' @param LDblocks.aux The LD block for the informative auxiliary population
#' @param LDblocks.tar The LD block for the target population
#' @param keep.test Participants to keep from the testing dataset
#' @param keep.ref Participants to keep from the reference panel
#'

prepare_data = function(ss_e, ss_a,pheno, ref.bfile.aux, test.bfile.aux,
                        ref.bfile.tar, sample_t,
                        test.bfile.tar, LDblocks.aux,
                        LDblocks.tar, keep.test = NULL,
                        keep.ref = NULL){
  cor_e <- lassosum::p2cor(p=ss_e$pvalue, n=ss_e$n, sign=ss_e$beta)

  adj = which(is.na(cor_e))
  if(length(adj)!=0){
    cor_e[adj] = sign(ss_e$beta[adj])*0.999
  }




  out_e <- lassosum::lassosum.pipeline(cor=cor_e,
                             chr=ss_e$chr,
                             pos=ss_e$bp,
                             A1=ss_e$effAllele,
                             A2=ss_e$refAllele,
                             ref.bfile = ref.bfile.aux,
                             test.bfile=test.bfile.aux,
                             keep.ref = keep.ref,
                             keep.test = keep.test,
                             LDblocks = LDblocks.aux)

  v_e = lassosum::validate(out_e, pheno = pheno)




  cor_a <- lassosum::p2cor(p=ss_a$pvalue, n=ss_a$n, sign=ss_a$beta)

  adj = which(is.na(cor_a))
  if(length(adj)!=0){
    cor_a[adj] = sign(ss_a$beta[adj])*0.999
  }


  out_a <- lassosum::lassosum.pipeline(cor=cor_a,
                             chr=ss_a$chr,
                             pos=ss_a$bp,
                             A1=ss_a$effAllele,
                             A2=ss_a$refAllele,
                             ref.bfile = ref.bfile.tar,
                             #keep.ref = sample_t,
                             keep.test = sample_t,
                             test.bfile = test.bfile.tar,
                             LDblocks = LDblocks.tar)
  v_a = lassosum::validate(out_a, pheno = pheno)

  return(list(validate.aux = v_e, validate.tar = v_a, output.aux = out_e,
              output.tar = out_a))
}


#' @title PRS_tf
#' @description  This function is used to apply TL-Multi to conduct multiethnic
#'   PRS prediction
#' @param ss_tl The GWAS data for TL-Multi
#' @param pheno The phenotype of the target population
#' @param ref.bfile The reference panel for the target population
#' @param test.bfile The testing panel for the target population
#' @param LDblocks The LD block for the target population
#' @param keep.test Participants to keep from the testing dataset
#' @param keep.ref Participants to keep from the reference panel
#'

PRS_tf = function(ss_tl, pheno, ref.bfile, test.bfile, LDblocks, keep.test = NULL,
                  keep.ref = NULL){

  cor_tl <- lassosum::p2cor(p=ss_tl$pvalue, n=ss_tl$n, sign=ss_tl$beta_tl)



  adj = which(is.na(cor_tl))
  if(length(adj)!=0){
    cor_tl[adj] = sign(ss_tl$beta_tl[adj])*0.999
  }




  out_tl <- lassosum::lassosum.pipeline(cor=cor_tl,
                              chr=ss_tl$chr,
                              pos=ss_tl$bp,
                              A1=ss_tl$effAllele,
                              A2=ss_tl$refAllele,
                              ref.bfile=ref.bfile,
                              test.bfile=test.bfile,
                              keep.ref = keep.ref,
                              keep.test = keep.test,
                              LDblocks = LDblocks)



  v_tl = lassosum::validate(out_tl, pheno=pheno)

  return(list(PRS = v_tl$best.pgs, best.beta = v_tl$best.beta))
}


#' @title ss_tl
#' @description  This function is used to generate the GWAS data for TL-Multi
#'   which were applied to conduct multiethnic PRS prediction
#' @details This function will return the multiethnic PRS calculated by TL-Multi
#' @param v_e The validation result of Lassosum for the auxiliary informative
#'   population
#' @param out_e The Lassosum result of the auxiliary informative population
#' @param ss_a The GWAS data for the target population
#' @param XTX The LD region of the target population
#' @param snp_list the rsid list of target population with colname 'snplist'
#' @return A list of multiethnic PRS calculated by TL-Multi.
#'
#'
#'


ss_tl <- function(v_e, out_e, ss_a, snp_list, XTX = XTX){


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

  ss_tl <-  merge(df_tl, snp_list, by.x="rsid", by.y="snplist", sort=F)
  ZTZ <-  XTX[ss_tl$rsid, ss_tl$rsid]

  ss_tl$beta_tl <-  ss_tl$beta_a - as.vector(ss_tl$beta_e %*% ZTZ)
  ss_tl$z <-  ss_tl$beta_tl / ss_tl$se
  ss_tl$pvalue <-  2 * stats::pnorm(-abs(ss_tl$z))





  return(ss_tl)

}


