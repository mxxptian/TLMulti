#' @title normalize
#' @description  This function is used to normalize genotypes matrix
#' @param x The genotype matrix of the target population
#' @return
#'
#'


normalize <- function(x){
  (x-mean(x)) / sd(x)
}
