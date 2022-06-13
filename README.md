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


result = PRS_tf(ss_e, ss_a, info_LD, tar_LD, Gt, PATH_TO_EUR, PATH_TO_REF, PATH_TO_TEST, pheno_t, sample_r = NULL,
                   sample_t = NULL, cluster = NULL)



```
