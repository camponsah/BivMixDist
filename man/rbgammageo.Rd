% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bgg_distribution.R
\name{rbgammageo}
\alias{rbgammageo}
\title{The bivariate gamma geometric distribution (BGG)}
\usage{
rbgammageo(n, alpha, beta, p)
}
\arguments{
\item{n}{size of sample.}

\item{alpha}{numeric parameter which must be greater than  0.}

\item{beta}{numeric parameter which must be  greater than  0.}

\item{p}{numeric parameter between 0 and 1.}
}
\value{
vector of random samples generate from BGG distribution.
}
\description{
Random sample generating function for bivariate gamma geometric distribution with parameters \eqn{\alpha >0, \beta > 0} and p in (0,1).
}
\details{
rbgammageo generates random samples from BGG distribution.
}
\examples{
data.df<-rbgammageo(20,alpha=1.5, beta=2, p=0.6)
data.df

}
\references{
Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Lévy process . Journal of Multivariate Analysis, 109:130-145.
\url{https://doi.org/10.1016/j.jmva.2012.03.004}
}
