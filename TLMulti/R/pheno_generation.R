#' @title pheno_generation
#' @description  This function is used to generate phenotypes of simulation data.
#' @details This function will return the phenotypes of two samples, in which most causal variants are shared. There is a difference between the target samples and the informative auxiliary samples in the genetic architecture, which causes estimation bias.
#' @param Za Genotypes of the training samples.
#' @param Ze Genotypes of the informative auxiliary samples.
#' @param Zt Genotypes of the testing samples.
#' @param ratio The proportion of causal SNPs to common SNPs.
#' @param rho The genetic correlation between these two populations.
#' @return A list of phenotypes according to the input genotypes.
#'
#' @examples
#' library(bigsnpr)
#' temp = snp_attachExtdata() # loading test data from 'bigsnpr' package
#' Ne = 300 # the number of samples for the informative auxiliary population used for simulating summary statistics
#' Na = 117 # the number of samples for the training target population used for simulating summary statistics
#' Ne = 100 # the individual-level data for the target population
#' Ge = temp$genotypes[1:300,]
#' Ga = temp$genotypes[301:417,]
#' Gt = temp$genotypes[418: 517,]
#' rho = 0.4 # the genetic architecture correlation
#' ratio = 0.01 # the causal snp proportion
#' pheno = pheno_generation(Ga, Ge, Gt, ratio, rho)




pheno_generation <- function(Za, Ze, Zt, ratio, rho){
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
