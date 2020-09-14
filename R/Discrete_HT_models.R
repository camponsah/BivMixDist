#' The family of discrete distributions with gamma hidden truncation.
#'
#' Random sample generating function for the family of discrete distributions with gamma hidden truncation: parameters \eqn{\alpha,\beta,\lambda,\delta, \mu \ge 0} and \eqn{p,\theta\in (0,1).} p in (0,1).
#'
#' rdHTgamma generates random sample from discrete Pareto distribution.
#'
#' @param n size of sample.
#' @param x  numeric truncation value which must be greater than or equal to  0.
#' @param alpha  numeric gamma shape parameter which must be greater than  0.
#' @param beta  numeric gamma rate parameter which must be greater than  0.
#' @param distr  a description of the discrete distribution which is truncated by gamma: "Poisson", "Geometric", "DPareto"(for discrete Pareto) and "nbinom" (for Negative binomial).
#' @param lambda numeric Poisson parameter which must be greater than 0.
#' @param delta  numeric discrete Pareto shape parameter which must be greater than or equal to 0.
#' @param prob numeric parameter between 0 and 1. Must be specified when distribution is either geometric or discrete Pareto or negative binomial.
#' @param size  numeric negative binomial parameter which must be greater than  0.
#'
#' @return  vector of samples generate from discrete Pareto distribution.
#'
#' @examples
#' N <- rdHTgamma(n=2, x=0.05, alpha = 1,beta = 100, distr = "Geometric",p = 0.5)
#' N
#'
#'@references  Amponsah, C. K. and  Kozubowski, T. J. (2021). Discrete Hidden Truncation models. In preparation.
#' \url{https://doi.org/not-done}
#'
#' @export
rdHTgamma<- function(n, x, alpha, beta, distr, lambda=NULL, prob=NULL,
                           delta=NULL, size=NULL){
  N <- NULL
  if(distr=="Poisson") {
    if (is.null(lambda)) stop("'lambda' must be specified")
    for (i in 1:n){
    K <- 0
    T <- 1
      while (K < T) {
      K <- stats:: rpois(n=1, lambda = lambda)
      gammaRv <- cumsum(stats:: rgamma(n=100000, shape = alpha, rate = beta))
      T <- sum(gammaRv <= x) + 1
      }
    N[i] <- K
    }
  }
  else if(distr=="Geometric") {
    if (is.null(prob)) stop("'prob' must be specified")
    for (i in 1:n){
      K <- 0
      T <- 1
      while (K < T) {
        K <- stats:: rgeom(n=1, prob = prob)
        gammaRv <- cumsum(stats:: rgamma(n=100000, shape = alpha, rate = beta))
        T <- sum(gammaRv <= x) + 1
      }
      N[i] <- K
    }
  }
  else if(distr=="DPareto") {
    if (is.null(delta) | is.null(prob)) stop("both 'delta' and 'prob' must be specified")
    for (i in 1:n){
      K <- 0
      T <- 1
      while (K < T) {
        K <- rdpareto(n=1, delta = delta, p=prob)
        gammaRv <- cumsum(stats:: rgamma(n=100000, shape = alpha, rate = beta))
        T <- sum(gammaRv <= x) + 1
      }
      N[i] <- K
    }
  }
  else if(distr=="nbinom") {
    if (is.null(size) | is.null(prob)) stop("both 'size' and 'prob' must be specified")
    for (i in 1:n){
      K <- 0
      T <- 1
      while (K < T) {
        K <- stats:: rnbinom(n=1,size =size , pron = prob)
        gammaRv <- cumsum(stats:: rgamma(n=100000, shape = alpha, rate = beta))
        T <- sum(gammaRv <= x) + 1
      }
      N[i] <- K
    }
  }
  return(N)
}



