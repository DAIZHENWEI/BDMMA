% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BDMMA.R
\name{fdr_cut}
\alias{fdr_cut}
\title{Threshold the posterior inclusion probability (PIP) through control Bayesian false discovery
rate (bFDR).}
\usage{
fdr_cut(PIP_vec, alpha = 0.1)
}
\arguments{
\item{PIP_vec}{A vector contains the PIPs of parameters}

\item{alpha}{The level of the bFDR to need to control (default = 0.1)}
}
\value{
The cutoff for PIPs to control the bFDR with the user defined value, alpha.
}
\description{
Threshold the posterior inclusion probability (PIP) through control Bayesian false discovery
rate (bFDR).
}
\examples{
data(L_mean)
cutoff <- fdr_cut(L_mean, alpha = 0.1)
}
