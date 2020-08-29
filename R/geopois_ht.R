#' The Geometric with Poisson hidden truncation distribution
#'
#' Random sample generating function for discrete Pareto distribution with parameters \eqn{\lambda > 0} and \eqn{p \in (0,1)}.
#'
#' rgeopois generates random sample from discrete geometric with Poisson hidden truncation distribution.
#'
#' @param n size of sample.
#' @param lambda  numeric parameter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of samples generate from discrete geometric with Poisson hidden truncation distribution.
#'
#' @examples
#' N <- rgeopois(20, lambda=2, p=0.3)
#' N
#'
#'@references  Amponsah, C. K. and  Kozubowski, T. J. (2020). Geometric with Poisson hidden truncation distribution.
#' \url{Yet to publish}
#'
#' @export
rgeopois <- function(n,lambda, p){
  K <- stats:: rgeom(n, prob = p)
  N <- stats:: rpois(n, lambda = lambda * (1-p))
  return(N+K)
}

#' The Geometric with Poisson hidden truncation distribution
#'
#' Probability mass function for geometric with Poisson hidden truncation distribution: parameters \eqn{\lambda > 0} and \eqn{p \in (0,1)}.
#'
#' dgeopois gives the probability mass.
#'
#' @param N A vector of random sample from discrete Pareto.
#' @param lambda  numeric parameter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return vector representing probabilities \eqn{P(N=n)}.
#'
#' @examples
#' prob<-dgeopois(seq(1,10,1), lambda=0.2, p=0.6)
#' prob
#'
#'@references  Amponsah, C. K. and  Kozubowski, T. J. (2020). Geometric with Poisson hidden truncation distribution.
#' \url{Yet to publish}
#'
#' @export
dgeopois <- function(N,lambda, p, log.p=FALSE){
  K <- stats:: ppois(N, lambda = lambda)
  W <-stats:: dgeom(N, prob = p) *exp(lambda *p)
  M <- K*W
  if (log.p == FALSE){
    return(M)
  }
  else {
    return(log(M))
  }
}

#' The Geometric with Poisson hidden truncation distribution
#'
#' Distribution function for Geometric with Poisson hidden truncation distribution, where parameters: \eqn{\lambda > 0} and \eqn{p \in (0,1)}.
#'
#' pgeopois gives the distribution function.
#'
#' @param q vector of quantiles representing numbers from discrete Pareto distribution.
#' @param lambda  numeric parameter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[N\leq n]}, otherwise, \eqn{$P[N> n]$}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return A vector cumulative probabilities
#'
#' @examples
#' prob<-pgeopois(seq(1,10, 1), lambda=0.2, p=0.6)
#' prob
#'
#'@references  Amponsah, C. K. and  Kozubowski, T. J. (2020). Geometric with Poisson hidden truncation distribution.
#' \url{Yet to publish}
#'
#' @export
pgeopois <- function(q,lambda, p, lower.tail=TRUE, log.p=FALSE){
  M <- NULL
  for (i in length(q)) {
    k <- seq(0, q[i])
    a <- 1- (1-p)^(q[i]-k+1)
    b <- (lambda *  (1-p))^(k)
    M[i]  <- exp(-lambda * (1-p)) * sum(a*b/gamma(k+1))
  }
  if (lower.tail == TRUE && log.p == FALSE){
    return(M)
  }
  else if (lower.tail == FALSE && log.p==FALSE){
    return(1-M)
  }
  else if (lower.tail == FALSE && log.p == TRUE){
    return(log(1-M))
  }
  else {
    return(log(M))
  }
}
