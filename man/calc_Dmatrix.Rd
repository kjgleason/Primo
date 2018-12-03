% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_Dmatrix}
\alias{calc_Dmatrix}
\title{Density matrix calculation}
\usage{
calc_Dmatrix(Q, D0, D1)
}
\arguments{
\item{Q}{matrix of configurations}

\item{D0}{estimate for the null density function}

\item{D1}{estimate for the alternative density function}
}
\value{
Returns a matrix estimating the conditional joint densities
under each configuration
(i.e. estimated probability of each configuration for each SNP).
}
\description{
Calculate conditional joint densities for each configuration, given marginal
null and alternative densities.
}