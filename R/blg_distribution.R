#' bivariate distribution with Lomax and geometric margins (BLG)
#'
#' Random sample generating function for bivariate distribution with Lomax and geometric margins (BLG) distribution with parameters \eqn{\alpha >0, \beta > 0} and p in [0,1].
#'
#' rblomaxgeo generates random samples from BLG distribution.
#'
#' @param n size of sample.
#' @param alpha paramter which must be numeric greater than  0.
#' @param beta  paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of samples generate from BLG distribution.
#'
#' @examples
#' data.df<-rblomaxgeo(20,alpha=1.5, beta=2, p=0.6)
#' data.df
#'
#'@references  Barreto-Souza, W. (2012). A bivariate distribution with Lomax and geometric margins . Journal of the Korean Statistical Society, 47:405-422.
#' \url{https://doi.org/10.1016/j.jkss.2018.04.006}
#'
#' @export
rblomaxgeo<- function(n,alpha,beta,p){
  N<- stats:: rgeom(n,p)+1
  X<- stats:: rgamma(n,shape =N,rate = beta )
  X<-X/stats::rgamma(1,shape =1/alpha,rate = 1/alpha )
  return(data.frame(X,N))
}

#' The bivariate distribution with Lomax and geometric margins (BLG)
#'
#' Density function for bivariate distribution with Lomax and geometric margins (BLG) distribution with parameters \eqn{\alpha >0, \beta > 0} and p in [0,1].
#'
#' dblomaxgeo is the density function for BLG model.
#'
#' @param data  bivariate vector  (X,N) observations from BLG model.
#' @param alpha paramter which must be numeric greater than  0.
#' @param beta  paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of densities.
#'
#' @examples
#' data.df<-rblomaxgeo(20, alpha=1.5, beta=2, p=0.6)
#' den<-dblomaxgeo(data.df, alpha=1.5, beta=2, p=0.6)
#' den
#'
#'@references  Barreto-Souza, W. (2012). A bivariate distribution with Lomax and geometric margins . Journal of the Korean Statistical Society, 47:405-422.
#' \url{https://doi.org/10.1016/j.jkss.2018.04.006}
#'
#' @export
dblomaxgeo<- function(data,alpha,beta,p,log.p=FALSE){
  N<-data[,2]
  X<-data[,1]
  M<-((1-p)^(N-1))*(p*beta^(alpha*N))*(X^(alpha*N-1))*exp(-beta*X)[3]
  if (log.p == FALSE){
    return(M)
  }else{
    M<- log(M)
    return(M)
  }
}


#' The bivariate distribution with Lomax and geometric margins (BLG)
#'
#' Density function for bivariate distribution with Lomax and geometric margins (BLG) distribution with parameters \eqn{\alpha >0, \beta > 0} and p in [0,1].
#'
#' pblomaxgeo is the distribution function for BLG model.
#'
#' @param data  bivariate vector  (X,N) observations from BLG model.
#' @param alpha paramter which must be numeric greater than  0.
#' @param beta  paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x, N \leq n]}, otherwise, \eqn{P[X > x, N > n]}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of distribution.
#'
#' @examples
#' data.df<-rblomaxgeo(20, alpha=1.5, beta=2, p=0.6)
#' den<-pblomaxgeo(data.df, alpha=1.5, beta=2, p=0.6)
#' den
#'
#'@references  Barreto-Souza, W. (2012). A bivariate distribution with Lomax and geometric margins . Journal of the Korean Statistical Society, 47:405-422.
#' \url{https://doi.org/10.1016/j.jkss.2018.04.006}
#'
#' @export
pblomaxgeo<- function(data,alpha,beta,p, lower.tail=TRUE,log.p=FALSE){
  N<-data[,2]
  X<-data[,1]
  M<-NULL
  t=1
  for (i in 1:length(N)) {
    k=seq(1,N[i])
    S0<-0
    for (j in k) {
      S0<-S0 + p*((1-p)^(j-1))*pracma::gammainc(beta*X[i],j*alpha)[3]
    }
    M[t]<- S0
    t=t+1
  }
  if (lower.tail==TRUE & log.p==FALSE){
    return(M)
  }
  else if (lower.tail==FALSE & log.p==TRUE){
    M<-log(1-M)
    return(M)
  }
  else if (lower.tail==TRUE & log.p==TRUE){
    M<-log(M)
    return(M)
  }
  else {
    return(1-M)
  }
}



