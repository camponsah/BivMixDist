#' The Discrete Pareto Distribution
#'
#' Random sample generation.
#'
#' rdpareto generates random sample from discrete pareto distribution.
#'
#' @param n size of sample.
#' @param delta  shape paramter which must be numeric greater than or equal to 0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of sample generate from discrete Pareto distribution.
#'
#' @examples
#' N<-rdpareto(20, delta=0.2,p=0.6)
#' N
#'
#'@references  Buddana, A. and  Kozubowski, T. J. (2014). Discrete Pareto distribution. Journal of Economics and Quality Control, 29(2):143-156.
#' \url{https://doi.org/10.1515/eqc-2014-0014}
#'
#' @export
rdpareto<-function(n,delta,p){
  u<- stats:: runif(n)
  sigma<- - 1/(delta*log(1 - p))
  return(ceiling(sigma*((1 - u)^(- delta) - 1)))
}

#' The Discrete Pareto Distribution
#'
#' Probaility mass function.
#'
#' ddpareto gives the probability mass.
#'
#' @param N A vector of random sample from discrete Pareto.
#' @param delta A shape paramter which must be numeric greater than or equal to 0.
#' @param p A numeric parameter between 0 and 1
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return vector representing probabilities \eqn{P(N=n)}.
#'
#' @examples
#' prob<-ddpareto(seq(1,10,1),delta=0.2,p=0.6)
#' prob
#'
#'@references  Buddana, A. and  Kozubowski, T. J. (2014). Discrete Pareto distribution. Journal of Economics and Quality Control, 29(2):143-156.
#' \url{https://doi.org/10.1515/eqc-2014-0014}
#'
#' @export
ddpareto<-function(N,delta,p,log.p=FALSE){
 W1<- (1 - delta*(N - 1)*log(1 - p))^( - 1/delta)
 W2<- (1 - delta*N*log(1 - p))^( - 1/delta)
 M<- W1 - W2
 if (log.p == FALSE){
   return(M)
   }else{
   M<- log(M)
   return(M)
   }
}

#' The Discrete Pareto Distribution
#'
#' Distribution function.
#'
#' pdpareto gives the distribution function.
#'
#' @param q vector of quantiles representing numbers from discrete pareto distribution.
#' @param delta  shape paramter which must be numeric greater than or equal to 0.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[N\leq n]}, otherwise, \eqn{$P[N> n]$}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return A vector cumulative probabilities
#'
#' @examples
#' prob<-pdpareto(seq(1,10,1),delta=0.2,p=0.6)
#' prob
#'
#'@references  Buddana, A. and  Kozubowski, T. J. (2014). Discrete Pareto distribution. Journal of Economics and Quality Control, 29(2):143-156.
#' \url{https://doi.org/10.1515/eqc-2014-0014}
#'
#' @export
pdpareto<- function(q,delta,p,lower.tail=TRUE,log.p=FALSE){
  sigma = -1/(delta*log(1 - p))
  M<- 1- ((1 + q/sigma)^(- 1/delta))
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

#' The Discrete Pareto Distribution
#'
#' Quantile function.
#'
#' qdpareto gives the quantile function.
#'
#' @param prob Vector of probabilities.
#' @param delta A shape paramter which must be numeric greater than or equal to 0.
#' @param p A numeric parameter between 0 and 1.
#'
#' @return  vector quantiles from discrete pareto distribution.
#'
#' @examples
#' q<-rdpareto(seq(0.1,0.6,0.1), delta=0.2,p=0.6)
#' q
#'
#'@references  Buddana, A. and  Kozubowski, T. J. (2014). Discrete Pareto distribution. Journal of Economics and Quality Control, 29(2):143-156.
#' \url{https://doi.org/10.1515/eqc-2014-0014}
#'
#' @export
qdpareto<-function(prob,delta,p){
  sigma <- - 1/(delta*log(1 - p))
  return(ceiling(sigma*((1-  prob)^(- delta) - 1)))
}

