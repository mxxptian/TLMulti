#' @title PRS_tf
#' @description  This function is used to apply TL-Multi to conduct multiethnic PRS prediction
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
                              ref.bfile=PATH_TO_DATA,
                              test.bfile=PATH_TO_DATA,
                              keep.ref = keep.ref,
                              keep.test = keep.test,
                              LDblocks = LDblocks)



  v_tl = validate(out_tl, pheno=tar_pheno)

  return(list(PRS = v_tl$best.pgs, best.beta = v_tl$best.beta))
}
