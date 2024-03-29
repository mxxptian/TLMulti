% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TLMulti.R
\name{prepare_data}
\alias{prepare_data}
\title{prepare_data}
\usage{
prepare_data(
  ss_e,
  ss_a,
  pheno,
  ref.bfile.aux,
  test.bfile.aux,
  ref.bfile.tar,
  sample_t,
  test.bfile.tar,
  LDblocks.aux,
  LDblocks.tar,
  keep.test = NULL,
  keep.ref = NULL
)
}
\arguments{
\item{ss_e}{The GWAS data for the informative auxiliary population}

\item{ss_a}{The GWAS data for the target population}

\item{pheno}{The phenotype of the target population}

\item{ref.bfile.aux}{The reference panel for the informative auxiliary
population}

\item{test.bfile.aux}{The testing panel for the informative auxiliary
population}

\item{ref.bfile.tar}{The reference panel for the target population}

\item{sample_t}{the test samples from data}

\item{test.bfile.tar}{The testing panel for the target population}

\item{LDblocks.aux}{The LD block for the informative auxiliary population}

\item{LDblocks.tar}{The LD block for the target population}

\item{keep.test}{Participants to keep from the testing dataset}

\item{keep.ref}{Participants to keep from the reference panel}
}
\description{
This function is used to prepare data for generating the GWAS
  data for TL-Multi. It was based on the lassosum proposed by Mak et al
  (2017).
}
