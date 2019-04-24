#' The bivariate gamma geometric distribution (BGG)
#'
#' Random sample generating function for bivariate gamma geometric distribution with parameters \eqn{\alpha >0, \beta > 0} and p in (0,1).
#'
#' rbgammageo generates random samples from BGG distribution.
#'
#' @param n size of sample.
#' @param alpha  shape paramter which must be numeric greater than  0.
#' @param beta  scale paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of samples generate from BGG distribution.
#'
#' @examples
#' data.df<-rbgammageo(20,alpha=1.5, beta=2, p=0.6)
#' data.df
#'
#'@references  Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Lévy process . Journal of Multivariate Analysis, 109:130-145.
#' \url{https://doi.org/10.1016/j.jmva.2012.03.004}
#'
#' @export
rbgammageo<- function(n,alpha,beta,p){
  N<- stats:: rgeom(n,p)+1
  X<- stats:: rgamma(n,shape =alpha*N,rate = beta )
  return(data.frame(X,N))
}


#' The bivariate gamma geometric distribution (BGG)
#'
#' Density function for bivariate gamma geometric distribution with parameters \eqn{\alpha >0, \beta > 0} and p in (0,1).
#'
#' dbgammageo is the density function for BGG model.
#'
#' @param data  bivariate vector  (X,N) observations from BGG model.
#' @param alpha  shape paramter which must be numeric greater than  0.
#' @param beta  scale paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of densities.
#'
#' @examples
#' data.df<-rbgammageo(20, alpha=1.5, beta=2, p=0.6)
#' den<-dbgammageo(data.df, alpha=1.5, beta=2, p=0.6)
#' den
#'
#'@references  Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Lévy process . Journal of Multivariate Analysis, 109:130-145.
#' \url{https://doi.org/10.1016/j.jmva.2012.03.004}
#'
#' @export
dbgammageo<- function(data,alpha,beta,p,log.p=FALSE){
  N<-data[,2]
  X<-data[,1]
  M<-((1-p)^(N-1))*(p*beta^(alpha*N))*(X^(alpha*N-1))*exp(-beta*X)/gamma(alpha*N)
  if (log.p == FALSE){
    return(M)
  }else{
    M<- log(M)
    return(M)
  }
}

#' The bivariate gamma geometric distribution (BGG)
#'
#' Distribution function for bivariate gamma geometric distribution with parameters \eqn{\alpha >0, \beta > 0} and p in (0,1).
#'
#' pbgammageo is the distribution function for BGG model.
#'
#' @param data  bivariate vector  (X,N) observations from BGG model.
#' @param alpha  shape paramter which must be numeric greater than  0.
#' @param beta  scale paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x, N \leq n]}, otherwise, \eqn{P[X > x, N > n]}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of distribution.
#'
#' @examples
#' data.df<-rbgammageo(20, alpha=1.5, beta=2, p=0.6)
#' den<-pbgammageo(data.df, alpha=1.5, beta=2, p=0.6)
#' den
#'
#'@references  Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Lévy process . Journal of Multivariate Analysis, 109:130-145.
#' \url{https://doi.org/10.1016/j.jmva.2012.03.004}
#'
#' @export
pbgammageo<- function(data,alpha,beta,p, lower.tail=TRUE,log.p=FALSE){
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

