# TLMulti

`TLMulti` is an R package to conduct multienthic polygenic risk scores (PRS). This algorithm borrows to the main ideas of transfer learning proposed in 2020 by  Li et al.(https://arxiv.org/abs/2006.10593) to extend Lassosum proposed by  Mak et al. in 2017 (https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.22050). The main challenge is that the majority of public genome-wide association study (GWAS) data has been conducted in European population. The accuracy of PRS prediction in non-European populations is diminished due to much smaller sample size of GWAS. TL-Multi treats the European population as informative auxiliary data and borrows the information to improve the learning performance of target population (e.g., non-European). TL-Multi only requests the summary statistics from European and the target populations and individual-level data from the target population. TL-Multi aims to improve the PRS prediction performance of the target population.

# Installation


You can install the development version of
`TLMulti` from Github via the `devtools` package. I suppose using
the `remotes` package would work as well.

``` r
devtools::install_github("mxxptian/TLMulti")
```

# Example

