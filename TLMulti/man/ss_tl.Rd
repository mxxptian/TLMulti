% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TLMulti.R
\name{ss_tl}
\alias{ss_tl}
\title{ss_tl}
\usage{
ss_tl(v_e, out_e, ss_a, snp_list, XTX = XTX)
}
\arguments{
\item{v_e}{The validation result of Lassosum for the auxiliary informative
population}

\item{out_e}{The Lassosum result of the auxiliary informative population}

\item{ss_a}{The GWAS data for the target population}

\item{snp_list}{the rsid list of target population with colname 'snplist'}

\item{XTX}{The LD region of the target population}
}
\value{
A list of multiethnic PRS calculated by TL-Multi.
}
\description{
This function is used to generate the GWAS data for TL-Multi
  which were applied to conduct multiethnic PRS prediction
}
\details{
This function will return the multiethnic PRS calculated by TL-Multi
}
