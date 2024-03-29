% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mult_mixed_distribution.R
\name{dmgammageo}
\alias{dmgammageo}
\title{The multivariate mixed gamma-geometric distributions (MMGG)}
\usage{
dmgammageo(data, alpha, beta, prob, theta, log.p = FALSE)
}
\arguments{
\item{data}{bivariate vector  (X,N) observations from BGG model.}

\item{alpha}{numeric vector of parameters which must be  greater than  0.}

\item{beta}{numeric vector of parameters which must be  greater than  0.}

\item{prob}{numeric vector parameters between 0 and 1.}

\item{theta}{numeric  parameter between -1 and (1+min(p))/(1-min(p)).}

\item{log.p}{logical; if TRUE, densities are given as log(den).}
}
\value{
vector of densities.
}
\description{
Density function for multivariate mixed gamma-geometric distribution with parameters \eqn{\alpha_1, \alpha_1 >0, \beta_1, \beta_1 > 0}, \eqn{p_1, P_2 \in (0,1)} and \eqn{-1\leq \theta \frac{1+\min(p_1,p_2)}{1+\min(p_1,p_2)}}.
}
\details{
dmgammageo is the density function.
}
\examples{
data.df<-rmgammageo(20, alpha=c(0.5, 1.5), beta=c(20, 50), prob=c(0.6, 0.4), theta= -0.2)
den<-dmgammageo(data.df, alpha=c(0.5, 1.5), beta=c(20, 50), prob=c(0.6, 0.4), theta= -0.2)
den

}
\references{
Amponsah, C. K. and Kozubowski, T.J., (2022). Inprint
}
