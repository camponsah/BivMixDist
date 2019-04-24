#' The bivariate exponential geometric distrubution (BEG)
#'
#' Random sample generating function for bivariate exponential geometric distribution with parameters \eqn{\beta > 0} and p in (0,1).
#'
#' rbexpgeo generates random sample from BEG distribution.
#'
#' @param n size of sample.
#' @param beta  scale paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of sample generate from BEG model.
#'
#' @examples
#' N<-rbexpgeo(20, beta=2,p=0.6)
#' N
#'
#'@references  Kozubowski, T.J., & Panorska, A.K. (2005). A Mixed bivariate distribution with exponential and geometric marginals. Journal of Statistical Planning and Inference, 134, 501-520.
#' \url{https://doi.org/10.1016/j.jspi.2004.04.010}
#'
#' @export
rbexpgeo<- function(n,beta,p){
  N<- stats:: rgeom(n,p)+1
  X<- stats:: rgamma(n,shape =N,rate = beta )
  return(data.frame(X,N))
}


#' The bivariate exponential geometric distrubution (BEG)
#'
#' Density function for bivariate exponential geometric distribution with parameters \eqn{\beta > 0} and p in (0,1).
#'
#' dbexpgeo is the density  function.
#'
#' @param data is bivariate vector  (X,N) vector representing observations from BEG model.
#' @param beta  scale paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of densities.
#'
#' @examples
#' data<-rbexpgeo(20, beta=2,p=0.6)
#' den<-dbexpgeo(data, beta=2,p=0.6)
#' den
#'
#'@references  Kozubowski, T.J., & Panorska, A.K. (2005). A Mixed bivariate distribution with exponential and geometric marginals. Journal of Statistical Planning and Inference, 134, 501-520.
#' \url{https://doi.org/10.1016/j.jspi.2004.04.010}
#'
#' @export
dbexpgeo<- function(data,beta,p,log.p=FALSE){
  N<-data[,2]
  X<-data[,1]
  M<-(p*beta^N)*((X*(1-p))^(N-1))*exp(-beta*X)/gamma(N)
  if (log.p == FALSE){
    return(M)
  }else{
    M<- log(M)
    return(M)
  }
}


#' The bivariate exponential geometric distrubution (BEG)
#'
#' Distribution function for bivariate exponential geometric distribution with parameters \eqn{\beta > 0} and p in (0,1).
#'
#' pbexpgeo is the distribution  function.
#'
#' @param data is bivariate vector  (X,N) vector representing observations from BEG model.
#' @param beta  scale paramter which must be numeric greater than 0.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x, N \leq n]}, otherwise, \eqn{P[X > x, N > n]}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of distribution.
#'
#' @examples
#' data<-rbexpgeo(20, beta=2,p=0.6)
#' den<-pbexpgeo(data, beta=2,p=0.6)
#' den
#'
#'@references  Kozubowski, T.J., & Panorska, A.K. (2005). A Mixed bivariate distribution with exponential and geometric marginals. Journal of Statistical Planning and Inference, 134, 501-520.
#' \url{https://doi.org/10.1016/j.jspi.2004.04.010}
#'
#' @export
pbexpgeo<- function(data,beta,p, lower.tail=TRUE,log.p=FALSE){
  N<-data[,2]
  X<-data[,1]
  S1<- stats:: ppois(N-1,lambda = ((1-p)*beta*X),lower.tail = TRUE)
  S2<- stats:: ppois(N-1,lambda = (beta*X),lower.tail = FALSE)
  M<-1-exp(-p*beta*X)*S1 -((1-p)^N)*S2
  if (lower.tail==FALSE & log.p==TRUE){
    M<-log(1-M)
    return(M)
  }
  else if (lower.tail==TRUE & log.p==TRUE){
    M<-log(M)
    return(M)
  }
  else if (lower.tail==TRUE & log.p==FALSE){
    return(M)
  }
  else {
    return(1-M)
  }
}


