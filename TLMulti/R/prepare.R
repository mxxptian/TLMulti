#' @title prepare_data
#' @description  This function is used to prepare data for generating the GWAS data for TL-Multi. It was based on the lassosum proposed by Mak et al (2017).
#' @param ss_e The GWAS data for the informative auxiliary population
#' @param ss_a The GWAS data for the target population
#' @param pheno The phenotype of the target population
#' @param ref.bfile.aux The reference panel for the informative auxiliary population
#' @param test.bfile.aux The testing panel for the informative auxiliary population
#' @param ref.bfile.tar The reference panel for the target population
#' @param test.bfile.tar The testing panel for the target population
#' @param LDblocks.aux The LD block for the informative auxiliary population
#' @param LDblocks.tar The LD block for the target population
#' @param keep.test Participants to keep from the testing dataset
#' @param keep.ref Participants to keep from the reference panel
#'

prepare_data = function(ss_e, ss_a,pheno, ref.bfile.aux, test.bfile.aux,
                        ref.bfile.tar,
                        test.bfile.tar, LDblocks.aux,
                        LDblocks.tar, keep.test = NULL,
                        keep.ref = NULL){
  cor_e <- p2cor(p=ss_e$pvalue, n=ss_e$n, sign=ss_e$beta)

  adj = which(is.na(cor_e))
  if(length(adj)!=0){
    cor_e[adj] = sign(ss_e$beta[adj])*0.999
  }




  out_e <- lassosum.pipeline(cor=cor_e,
                             chr=ss_e$chr,
                             pos=ss_e$bp,
                             A1=ss_e$effAllele,
                             A2=ss_e$refAllele,
                             ref.bfile = ref.bfile.aux,
                             test.bfile=test.bfile.aux,
                             keep.ref = keep.ref,
                             keep.test = keep.test,
                             LDblocks = LDblocks.aux)

  v_e = validate(out_e, pheno=tar_pheno)

  tmp_dir = tempdir()
  file_name = list.files(tmp_dir)

  if(length(file_name[grep('.bk', file_name)])>=1){
    file.remove(paste0(tmp_dir,'/', file_name[grep('.bk', file_name)]))

  }




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
                             ref.bfile=PATH_TO_DATA,
                             #keep.ref = sample_t,
                             keep.test = sample_t,
                             test.bfile=PATH_TO_DATA,
                             LDblocks = tar_LD)
  v_a = validate(out_a, pheno=tar_pheno)

  return(list(validate.aux = v_e, validate.tar = v_a, output.aux = out_e,
              output.tar = out_a))
}
