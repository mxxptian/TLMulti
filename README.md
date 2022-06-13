# TLMulti

`TLMulti` is an R package to conduct multienthic polygenic risk scores (PRS). This algorithm borrows to the main ideas of transfer learning proposed in 2020 by  Li et al.(https://arxiv.org/abs/2006.10593) to extend Lassosum proposed by  Mak et al. in 2017 (https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.22050). The main challenge is that the majority of public genome-wide association study (GWAS) data has been conducted in European population. The accuracy of PRS prediction in non-European populations is diminished due to much smaller sample size of GWAS. TL-Multi treats the European population as informative auxiliary data and borrows the information to improve the learning performance of target population (e.g., non-European). TL-Multi only requests the summary statistics from European and the target populations and individual-level data from the target population. TL-Multi aims to improve the PRS prediction performance of the target population.

# Installation


You can install the development version of
`TLMulti` from Github via the `devtools` package. I suppose using
the `remotes` package would work as well.

Before installation of TL-Multi, you are also requested the below packages:
``` r
install.packages(c('bigsnpr', 'bigstatsr', 'lassosum', 'data.table', 'parallel', 'dplyr', 'mvtnorm'), dependencies=TRUE)

```

``` r
devtools::install_github("mxxptian/TLMulti")
```

# Example

``` r

library(bigsnpr)
library(mvtnorm)
library(lassosum)
library(genio)

temp = snp_attachExtdata()
map = temp$map

Ne = 300 
Na = 117 
Nt = 100 

Ge = temp$genotypes[1:300,]
Ga = temp$genotypes[301:417,]
Gt = temp$genotypes[418: 517,] 


Ze <- apply(Ge, 2, normalize)  # normalized genotype for european
Za <- apply(Ga, 2, normalize)  # normalized genotype for asian
Zt <- apply(Gt, 2, normalize)  # normalized genotype for test

rho = 0.4
h2 = 0.5
ratio = 0.01
pheno = pheno_generation(Ne, Na, Nt, Za, Ze, Zt, ratio, rho, h2)


pheno_a <- pheno$pheno_a
pheno_e <- pheno$pheno_e
pheno_t <- pheno$pheno_t

estim_e <- big_univLinReg(X = as_FBM(Ze), pheno_e)  # GWAS for EUR
ss_e <- cbind(map, estim_e)
names(ss_e) = c("chr", "rsid", "genetic.dist", "bp", "effAllele", "refAllele",
                "beta", "se", "z")
ss_e$pvalue = 2*pnorm(-abs(ss_e$z))
ss_e$n <- Ne


estim_a <- big_univLinReg(X = as_FBM(Za), pheno_a)  # GWAS for ASN
ss_a <- cbind(map, estim_a)
names(ss_a) = c("chr", "rsid", "genetic.dist", "bp", "effAllele", "refAllele",
                "beta", "se", "z")
ss_a$pvalue = 2*pnorm(-abs(ss_a$z))
ss_a$n <- Na

PATH_TO_TEST =  'testsample' #The test panel corresponding to the target population (PLINK files)
PATH_TO_REF = "1000G_EAS.chr1-22" #the reference panel for the target population (PLINK file)
PATH_TO_EUR = '1000G_EUR_Phase3_hm3only1.2all' #The reference panel for informative auxiliary population (PLINK file)
info_LD = 'EUR.hg38' #LD region file for the informative auxiliary population
tar_LD = 'ASN.hg38' #LD region file for the target population


result = PRS_tf(ss_e, ss_a, info_LD, tar_LD, X,PATH_TO_EUR, PATH_TO_REF, PATH_TO_TEST, tar_pheno, sample_r = NULL,
                   sample_t = NULL){
  
  Z = apply(Gt, 2, normalize)
  XTX = t(Z)%*%Z/nrow(Z)
  
  
  cor_e = p2cor(p=ss_e$pvalue, n=ss_e$n, sign=ss_e$beta)
  
  
  out_e <- lassosum.pipeline(cor=cor_e,
                             chr=ss_e$chr,
                             pos=ss_e$bp,
                             A1=ss_e$effAllele,
                             A2=ss_e$refAllele,
                             ref.bfile = PATH_TO_EUR,
                             test.bfile=PATH_TO_TEST,
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





```
