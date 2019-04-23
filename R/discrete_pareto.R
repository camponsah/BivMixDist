#' The Discrete Pareto Distribution
#'
#' Random sample generation.
#' @usage
#' @param n Size of sample.
#' @param delta A shape paramter which must be numeric greater than or equal to 0.
#' @param p A numeric parameter between 0 and 1
#' @return A vector of sample generate from discrete Pareto distribution.
#' @export
rdpareto<-function(n,delta,p){
  u<- runif(n)
  sigma<- - 1/(delta*log(1 - p))
  return(ceiling(sigma*((1 - u)^(- delta) - 1)))
}

#' The Discrete Pareto Distribution
#'
#' Probaility mass function.
#' @usage
#' @param N A vector of random sample from discrete Pareto.
#' @param delta A shape paramter which must be numeric greater than or equal to 0.
#' @param p A numeric parameter between 0 and 1
#' @return A vector consiting of densities
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
#' @usage
#' @param q A vector of quantiles values from discrete pareto distribution.
#' @param delta A shape paramter which must be numeric greater than or equal to 0.
#' @param p A numeric parameter between 0 and 1
#' @param lower.tail cummulative probabilities if TRUE otherwise survival probilities
#' @param log.p log probabilities if TRUE
#' @return A vector cumulative probabilities
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
#' @usage
#' @param prob Vector of probabilities.
#' @param delta A shape paramter which must be numeric greater than or equal to 0.
#' @param p A numeric parameter between 0 and 1
#' @return A vector quantiles
#' @export
qdpareto<-function(prob,delta,p){
  sigma <- - 1/(delta*log(1 - p))
  return(ceiling(sigma*((1-  prob)^(- delta) - 1)))
}

