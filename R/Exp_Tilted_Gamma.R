
#' The Exponentially Tilted Gamma (ETG) distribution
#'
#' Probability density function of ETG with parameters \eqn{\alpha,\beta \ge 0} and \eqn{\kappa \in [0,1]}.
#'
#' detgamma gives the probability density(ies).
#'
#' @param X A vector of random sample from ETG.
#' @param alpha A shape parameter which must be numeric greater than or equal to 0.
#' @param beta A scale parameter which must be numeric greater than or equal to 0.
#' @param K A numeric parameter between 0 and 1
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return vector representing densities \eqn{P(X=x)}.
#'
#' @examples
#' prob <- detgamma(seq(0.1,5,0.1),alpha=0.2,beta=10, K=0.2)
#' prob
#'
#'@references  Amponsah, C. K. and  Kozubowski, T. J. (2021). A computational approach to estimation of discrete Pareto parameters. Preprint
#'
#' @export
detgamma <- function(X, alpha, beta, K,log.p=FALSE){
  w <- K^alpha
  M1 <- stats:: dgamma(X, shape = alpha, rate = beta) * (1/(1-w))
  M2 <- stats:: dgamma(X, shape = alpha, rate = beta/K) * (w/(1-w))
  M <- M1 - M2
  if (log.p == FALSE){
    return(M)
  }else{
    M<- log(M)
    return(M)
  }
}

#' The Exponentially Tilted Gamma (ETG) distribution
#'
#' Probability distribution function of ETG with parameters \eqn{\alpha,\beta \ge 0} and \eqn{\kappa \in [0,1]}.
#'
#' petgamma gives the probability distribution.
#'
#' @param X A vector of random sample from ETG.
#' @param alpha A shape parameter which must be numeric greater than or equal to 0.
#' @param beta A scale parameter which must be numeric greater than or equal to 0.
#' @param K A numeric parameter between 0 and 1
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X\leq n)}, otherwise, \eqn{$P[N> n]$}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return vector representing densities \eqn{P(X=x)}.
#'
#' @examples
#' prob <- petgamma(seq(0.1,1,0.1),alpha=0.2,beta=10, K=0.2)
#' prob
#'
#'@references  Amponsah, C. K. and  Kozubowski, T. J. (2021). A computational approach to estimation of discrete Pareto parameters. Preprint
#'
#' @export
petgamma <- function(X, alpha, beta, K,lower.tail=TRUE,log.p=FALSE){
  w <- K^alpha
  M1 <- stats:: pgamma(X, shape = alpha, rate = beta) * (1/(1-w))
  M2 <- stats:: pgamma(X, shape = alpha, rate = beta/K) * (w/(1-w))
  M <- M1 - M2
  if (lower.tail == TRUE & log.p == FALSE){
    return(M)
  }else if (lower.tail == TRUE & log.p == TRUE){
    M<- log(M)
    return(M)
  }else if (lower.tail == FALSE & log.p == TRUE){
    M<- log(1 - M)
    return(M)
  }else {
    return(1 - M)
  }
}
