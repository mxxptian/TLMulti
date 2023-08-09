# TLMulti

`TLMulti` is an R package to conduct multienthic polygenic risk scores (PRS). This algorithm borrows to the main ideas of transfer learning proposed in 2020 by  Li et al.(https://arxiv.org/abs/2006.10593) to extend Lassosum proposed by  Mak et al. in 2017 (https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.22050). The main challenge is that the majority of public genome-wide association study (GWAS) data has been conducted in European population. The accuracy of PRS prediction in non-European populations is diminished due to much smaller sample size of GWAS. TL-Multi treats the European population as informative auxiliary data and borrows the information to improve the learning performance of target population (e.g., non-European). TL-Multi only requests the summary statistics from European and the target populations and individual-level data from the target population. TL-Multi aims to improve the PRS prediction performance of the target population.

# Installation


You can install the development version of
`TLMulti` from Github via the `devtools` package. I suppose using
the `remotes` package would work as well.

Before installation of TL-Multi, you are also requested the below packages:
``` r
install.packages(c('lassosum', 'mvtnorm'), dependencies=TRUE)

```

``` r
devtools::install_github("mxxptian/TLMulti")
```
Or you can also install by the source file:

``` r
install.packages("path/TLMulti_0.1.0.tar.gz", repos=NULL)
```


# Example

``` r
#define a function to simulate phenotype

pheno_generation <- function(Ne, Na, Nt, Za, Ze, Zt, ratio, rho, h2){

  M <- dim(Za)[2]  # number of SNPs
  m <- ceiling(M*ratio)  # number of causal SNPs

  set <- sample(1:M, m)  # index for causal SNPs

  beta_e <- rep(0, M)  # coef for european
  beta_a <- rep(0, M)  # coef for asian


  b <- rmvnorm(m, sigma = matrix(data = h2/m*c(1, rho, rho, 1), nrow = 2))
  beta_e[set] <- b[,1]
  beta_a[set] <- b[,2]

  pheno_e <- as.vector(Ze%*%beta_e+rnorm(Ne, 0, sqrt(1-h2)))  # phenotype for european
  pheno_a <- as.vector(Za%*%beta_a+rnorm(Na, 0, sqrt(1-h2)))  # phenotype for asian
  pheno_t <- as.vector(Zt%*%beta_a+rnorm(Nt, 0, sqrt(1-h2)))  # phenotype for test

  return(list(pheno_a=pheno_a, pheno_e=pheno_e, pheno_t=pheno_t))
}


library(bigsnpr)
library(mvtnorm)
library(lassosum)
library(genio)

G = snp_attach('./HK_impute/HK.chr1-22.impute_QC.rds')
Geur <- snp_attach("ukb_imp/ukbbk.rds")
temp = table(Geur$map$rsid)
temp2 = names(temp)[temp>=2]
rm(temp)


Geur$genotypes = Geur$genotypes[Ne,!(Geur$map$rsid %in% temp2)]
Geur$map = Geur$map[!(Geur$map$rsid %in% temp2),]

common_snp = intersect(Geur$map$rsid, G$map$marker.ID)


map_hk = G$map[(G$map$chromosome %in%chr_asn) & (G$map$marker.ID %in% common_snp),]


ga = G$genotypes[,(G$map$chromosome %in%chr_asn) & (G$map$marker.ID %in% common_snp)]
ge = Geur$genotypes[,(Geur$map$chromosome%in%chr_eur)& (Geur$map$rsid %in% common_snp)]


Ge = ge
Ga = ga[1:Na,]
Gt = ga[(Na+1):(Na+Nt),]



Ze <- apply(Ge, 2, normalize)  # normalized genotype for european
Za <- apply(Ga, 2, normalize)  # normalized genotype for asian
Zt <- apply(Gt, 2, normalize)  # normalized genotype for test

map <- map_hk

sample_t <- rep(FALSE, nrow(G$genotypes)); sample_t[(Na+1):(Na+Nt)] <- TRUE  # sample of test



pheno = pheno_generation(Ne, Na, Nt, Za, Ze, Zt,
                          ratio = 0.01, rho = 0.2, h2 = 0.5)


pheno_a <- pheno$pheno_a
pheno_e <- pheno$pheno_e
pheno_t <- pheno$pheno_t

XTX <- t(Zt)%*%Zt/Nt
snplist <-map$marker.ID
rownames(XTX) = snplist
colnames(XTX) = snplist
snplist =  data.frame(snplist)

ss_e <- big_univLinReg(X = as_FBM(Ze), pheno_e)  # GWAS for EUR
ss_e <- cbind(map, ss_e)
names(ss_e) = c("chr", "rsid", "genetic.dist", "bp", "effAllele", "refAllele",
                "beta", "se", "z")
ss_e$pvalue = 2*pnorm(-abs(ss_e$z))
ss_e$n <- Ne

ss_a <- big_univLinReg(X = as_FBM(Za), pheno_a)  # GWAS for ASN
ss_a <- cbind(map, ss_a)
names(ss_a) = c("chr", "rsid", "genetic.dist", "bp", "effAllele", "refAllele",
                "beta", "se", "z")
ss_a$pvalue = 2*pnorm(-abs(ss_a$z))
ss_a$n <- Na


info_LD = 'EUR.hg38'
tar_LD = 'ASN.hg38'
PATH_TO_DATA = './HK_impute/HK.chr1-22.impute'

sample_r = NULL
sample_t = sample_t
cluster = NULL


pheno_t_table <- G$fam[sample_t,c(1,2)]


colnames(pheno_t_table) <- c('FID','IID')
pheno_t_table$pheno <- pheno_t # FAM for test

tar_pheno = pheno_t_table


pre = prepare_data(ss_e, ss_a, pheno = tar_pheno, ref.bfile.aux = PATH_TO_DATA, test.bfile.aux =PATH_TO_DATA,
                        ref.bfile.tar =PATH_TO_DATA,
                        test.bfile.tar=PATH_TO_DATA, LDblocks.aux = 'EUR.hg38',
                        LDblocks.tar = 'ASN.hg38', keep.test = sample_t,
                        keep.ref = NULL)

ss_tl = ss_tl(v_e = pre$validate.aux, out_e = pre$output.aux,
              ss_a = ss_a, XTX = XTX)

result = PRS_tf(ss_tl, pheno = tar_pheno, ref.bfile = PATH_TO_DATA, test.bfile = PATH_TO_DATA,
                                LDblocks  = 'ASN.hg38', keep.ref = NULL, keep.test = sample_t)




```
# Citation
Tian, P., Chan, T. H., Wang, Y.-F., Yang, W., Yin, G., and Zhang, Y. D. (2022). [Multiethnic polygenic risk prediction in diverse populations through transfer learning.](http://journal.frontiersin.org/article/10.3389/fgene.2022.906965/full?&utm_source=Email_to_authors_&utm_medium=Email&utm_content=T1_11.5e1_author&utm_campaign=Email_publication&field=&journalName=Frontiers_in_Genetics&id=906965) Frontiers in Genetics, 13.
