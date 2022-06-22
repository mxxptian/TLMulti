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

G = snp_attach('./HK_impute/HK.chr1-22.impute_QC.rds') # loading the genotype for the target population
Geur <- snp_attach("ukb_imp/ukbbk.rds") # loading the genotype for the auxiliary informative population 
temp = table(Geur$map$rsid)
temp2 = names(temp)[temp>=2]
rm(temp)


Geur$genotypes = Geur$genotypes[Ne,!(Geur$map$rsid %in% temp2)]
Geur$map = Geur$map[!(Geur$map$rsid %in% temp2),] # remove the duplicated SNPs

common_snp = intersect(Geur$map$rsid, G$map$marker.ID) # draw the common SNPs between the target population and the auxiliary informative population


map_hk = G$map[(G$map$chromosome %in%chr_asn) & (G$map$marker.ID %in% common_snp),]


ga = G$genotypes[,(G$map$chromosome %in%chr_asn) & (G$map$marker.ID %in% common_snp)] # sample the common SNPs for the target population
ge = Geur$genotypes[,(Geur$map$chromosome%in%chr_eur)& (Geur$map$rsid %in% common_snp)] # sample the common SNPs for the auxiliary informative population


Ge = ge # training data of the auxiliary informative population
Ga = ga[1:Na,] # training data of the target population
Gt = ga[(Na+1):(Na+Nt),] # testing data of the target population



Ze <- apply(Ge, 2, normalize)  # normalized genotype for european
Za <- apply(Ga, 2, normalize)  # normalized genotype for asian
Zt <- apply(Gt, 2, normalize)  # normalized genotype for test

map <- map_hk

sample_t <- rep(FALSE, nrow(G$genotypes)); sample_t[(Na+1):(Na+Nt)] <- TRUE  # sample of test



pheno = pheno_generation(Ne, Na, Nt, Za, Ze, Zt, 
                          ratio = 0.01, rho = 0.2, h2 = 0.5) #simulate the corresponding phenotype


pheno_a <- pheno$pheno_a # phenotype of the training data for the target population
pheno_e <- pheno$pheno_e # phenotype for the auxiliary informative population
pheno_t <- pheno$pheno_t # phenotype of the testing data for the target population 

XTX <- t(Zt)%*%Zt/Nt # compute the LD of the target population
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

tar_pheno = pheno_t_table #generate testing phenotype data


out_e <- lassosum.pipeline(cor=cor_e, 
                           chr=ss_e$chr, 
                           pos=ss_e$bp, 
                           A1=ss_e$effAllele, 
                           A2=ss_e$refAllele, 
                           ref.bfile = PATH_TO_DATA, 
                           test.bfile=PATH_TO_DATA,
                           #keep.ref = sample_t,
                           keep.test = sample_t,
                           LDblocks = info_LD)

v_e = validate(out_e, pheno=tar_pheno)

t12 = Sys.time()

cor_a <- p2cor(p=ss_a$pvalue, n=ss_a$n, sign=ss_a$beta)


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


ss_tl = ss_tl(v_e = v_e, out_e = out_e, ss_a = ss_a, XTX = XTX) # generate the GWAS for TL-Multi

cor_tl <- p2cor(p=ss_tl$pvalue, n=ss_tl$n, sign=ss_tl$beta_tl)


out_tl <- lassosum.pipeline(cor=cor_tl, 
                            chr=ss_tl$chr, 
                            pos=ss_tl$bp, 
                            A1=ss_tl$effAllele, 
                            A2=ss_tl$refAllele, 
                            ref.bfile=PATH_TO_DATA, 
                            test.bfile=PATH_TO_DATA, 
                            keep.test = sample_t,
                            LDblocks = 'ASN.hg38')



v_tl = validate(out_tl, pheno=tar_pheno)



```
