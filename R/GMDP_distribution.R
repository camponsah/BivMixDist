#' The gamma mixture distcrete Pareto distribution (GMDP)
#'
#' Random sample generating function for gamma mixture distcrete Pareto distribution with parameters \eqn{\alpha >0, \beta > 0, \delta \ge 0} and p in (0,1).
#'
#' rgammamixdpareto generates random samples from GMDP distribution.
#'
#' @param n size of sample.
#' @param alpha  numeric paramter which must be greater than  0.
#' @param beta  numeric paramter which must be greater than  0.
#' @param  delta numeric paramter which must be greater than or equal to  0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of random samples from GMDP distribution.
#'
#' @examples
#' data.df<-rgammamixdpareto(20,alpha=1.5, beta=2, delta=0.1, p=0.3)
#' data.df
#'
#'@references  Amponsah, C. K., Kozubowski, T. J. and Panorska, A. K. (2019). A bivariate distribution with gamma mixture discrete Pareto marginals . In print.
#' \url{https://scholarworks.unr.edu/bitstream/handle/11714/2065/Amponsah_unr_0139M_12378.pdf?sequence=1&isAllowed=y}
#'
#' @export
rgammamixdpareto<- function(n,alpha,beta,delta,p){
  u<- stats:: runif(n)
  sigma=-1/(delta*log(1-p))
  N<-ceiling(sigma*((1-u)^(-delta) -1))
  X<- stats:: rgamma(n,shape =alpha*N,rate = beta )
  return(data.frame(X,N))
}

#' The gamma mixture distcrete Pareto distribution (GMDP)
#'
#' Density function for gamma mixture distcrete Pareto distribution with parameters \eqn{\alpha >0, \beta > 0, \delta \ge 0} and p in (0,1).
#'
#' dgammamixdpareto is the density function of GMDP model.
#'
#' @param data  vector  (X,N) representing observation from GMDP model.
#' @param alpha  numeric paramter which must be greater than  0.
#' @param beta  numeric paramter which must be greater than  0.
#' @param  delta numeric paramter which must be greater than or equal to  0.
#' @param p numeric parameter between 0 and 1.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of densities.
#'
#' @examples
#' data.df<-rgammamixdpareto(20,alpha=1.5, beta=2, delta=0.1, p=0.3)
#' den<- dgammamixdpareto(data.df,alpha=1.5, beta=2, delta=0.1, p=0.3)
#' den
#'
#'@references  Amponsah, C. K., Kozubowski, T. J. and Panorska, A. K. (2019). A bivariate distribution with gamma mixture discrete Pareto marginals . In print.
#' \url{https://scholarworks.unr.edu/bitstream/handle/11714/2065/Amponsah_unr_0139M_12378.pdf?sequence=1&isAllowed=y}
#'
#' @export
dgammamixdpareto<- function(data,alpha,beta,delta,p,log.p=FALSE){
  N<-data[,2]
  X<-data[,1]
  p1<-(beta^(alpha*N))*(X^(alpha*N-1))*exp(-beta*X)/gamma(alpha*N)
  p21<-(1-delta*(N-1)*log(1-p))^(-1/delta)
  p22<-(1-delta*N*log(1-p))^(-1/delta)
  M<-p1*(p21-p22)
  if (log.p == FALSE){
    return(M)
  }else{
    M<- log(M)
    return(M)
  }
}

#' The gamma mixture distcrete Pareto distribution (GMDP)
#'
#' Distribution function for gamma mixture distcrete Pareto distribution with parameters \eqn{\alpha >0, \beta > 0, \delta \ge 0} and p in (0,1).
#'
#' pgammamixdpareto is the distribution function of GMDP model.
#'
#' @param data  vector  (X,N) representing observation from GMDP model.
#' @param alpha  numeric paramter which must be greater than  0.
#' @param beta  numeric paramter which must be greater than  0.
#' @param  delta numeric paramter which must be greater than or equal to  0.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x, N \leq n]}, otherwise, \eqn{P[X > x, N > n]}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of probabilities.
#'
#' @examples
#' data.df<-rgammamixdpareto(20,alpha=1.5, beta=2, delta=0.1, p=0.3)
#' prob<- pgammamixdpareto(data.df,alpha=1.5, beta=2, delta=0.1, p=0.3)
#' prob
#'
#'@references  Amponsah, C. K., Kozubowski, T. J. and Panorska, A. K. (2019). A bivariate distribution with gamma mixture discrete Pareto marginals . In print.
#' \url{https://scholarworks.unr.edu/bitstream/handle/11714/2065/Amponsah_unr_0139M_12378.pdf?sequence=1&isAllowed=y}
#'
#' @export
pgammamixdpareto<- function(data,alpha,beta,delta,p, lower.tail=TRUE,log.p=FALSE){
  N<-data[,2]
  X<-data[,1]
  M<-NULL
  t=1
  for (i in 1:length(N)) {
    k=seq(1,N[i])
    S0<-0
    for (j in k) {
      p21<-(1-delta*(N[i]-1)*log(1-p))^(-1/delta)
      p22<-(1-delta*N[i]*log(1-p))^(-1/delta)
      S0<-S0 + (p21-p22)*pracma::gammainc(beta*X[i],j*alpha)[3]
    }
    M[t]<- S0
    t=t+1
  }
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

