% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixbeg_EM.R
\name{func_mgeo_em}
\alias{func_mgeo_em}
\title{E and M steps for mxiture of geometric part of FMBEG}
\usage{
func_mgeo_em(data, p, pi, l)
}
\arguments{
\item{data}{data-vector from FMBEG distribution.}

\item{p}{vector of numeric parameters between 0 and 1.}

\item{pi}{vector of geometric membership probabilities between 0 and 1, and sum equal to 1.}

\item{l}{number of geometric components.}
}
\value{
list containing  parameter estimates of p and pi
}
\description{
This function computes the p and pi of FMBEG distribution.
}
