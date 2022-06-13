#' @title PRS_tf
#' @description  This function is used to run TL-Multi which is applied to conduct multiethnic PRS prediction
#' @details This function will return the multiethnic PRS calculated by TL-Multi
#' @param ss_a GWAS of the target samples.
#' @param ss_e GWAS of the informative auxiliary samples.
#' @param info_LD LD matrix in Lassosum for the informative auxiliary samples
#' @param tar_LD LD matrix in Lassosum for the target samples
#' @param X The individual-level data of the target samples.
#' @param PATH_TO_REF The data path of reference panel and test panel in Lassosum.
#' @param PATH_TO_TEST The data path of reference panel and test panel in Lassosum.
#' @param sample_r Participants to keep from the reference panel
#' @param sample_t Participants to keep from the test panel
#' @return A list of multiethnic PRS calculated by TL-Multi.
#'
#'
PRS_tf <- function(cor_e, cor_a, info_LD, tar_LD, X, PATH_TO_REF, PATH_TO_TEST, tar_pheno, sample_r = NULL,
                   sample_t = NULL){

  Z = apply(X, 2, normalize)
  XTX = t(Z)%*%Z/nrow(Z)





  out_e <- lassosum.pipeline(cor=cor_e,
                             chr=ss_e$chr,
                             pos=ss_e$bp,
                             A1=ss_e$effAllele,
                             A2=ss_e$refAllele,
                             ref.bfile = PATH_TO_TEST,
                             test.bfile=PATH_TO_REF,
                             keep.ref = sample_r,
                             keep.test = sample_t,
                             LDblocks = info_LD)
  v_e = validate(out_e, tar_pheno)

  t12 = Sys.time()

  cor_a <- p2cor(p=ss_a$pvalue, n=ss_a$n, sign=ss_a$beta)

  adj = which(is.na(cor_a))
  if(length(adj)!=0){
    cor_a[adj] = sign(ss_a$beta[adj])*0.999
  }



  out_a <- lassosum.pipeline(cor=cor_a,
                             chr=ss_a$chr,
                             pos=ss_a$bp,
                             A1=ss_a$effAllele,
                             A2=ss_a$refAllele,
                             test.bfile = PATH_TO_TEST,
                             ref.bfile=PATH_TO_REF,
                             keep.ref = sample_r,
                             keep.test = sample_t,
                             test.bfile=PATH_TO_DATA,
                             LDblocks = tar_LD)

  v_a = validate(out_a, pheno=tar_pheno)

  ##### Conduct Lassosum for Transfer Learning


  best_beta <- data.frame(v_e$best.beta, out_e$sumstats$pos, out_e$sumstats$chr)
  colnames(best_beta) <- c('beta_e', 'bp', 'chr')
  df_tl <-  ss_a
  names(df_tl)[names(df_tl) == "beta"] = "beta_a"
  df_tl <-  merge(df_tl, best_beta[index,], by.x=c('bp','chr'), by.y=c("bp","chr"),  sort=F)

  ss_tl <-  merge(df_tl, snplist, by.x="rsid", by.y="snplist", sort=F)
  ZTZ <-  XTX[ss_tl$rsid, ss_tl$rsid]

  ss_tl$beta_tl <-  ss_tl$beta_a - as.vector(ss_tl$beta_e %*% ZTZ)
  ss_tl$z <-  ss_tl$beta_tl / ss_tl$se
  ss_tl$pvalue <-  2 * pnorm(-abs(ss_tl$z))


  cor_tl <- p2cor(p=ss_tl$pvalue, n=ss_tl$n, sign=ss_tl$beta_tl)


  adj = which(is.na(cor_tl))
  if(length(adj)!=0){
    cor_tl[adj] = sign(ss_tl$beta_tl[adj])*0.999
  }

  out_tl <- lassosum.pipeline(cor=cor_tl,
                              chr=ss_tl$chr,
                              pos=ss_tl$bp,
                              A1=ss_tl$effAllele,
                              A2=ss_tl$refAllele,
                              ref.bfile=PATH_TO_REF,
                              test.bfile=PATH_TO_TEST,
                              keep.ref = sample_r,
                              keep.test = sample_t,
                              LDblocks = tar_LD)


  v_tl = validate(out_tl, pheno=tar_pheno)


  R2_tl <- cor(v_e$best.pgs+v_tl$best.pgs,pheno_t)^2


  return(list(PRS_tl=v_e$best.pgs+v_tl$best.pgs))

}


