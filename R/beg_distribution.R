#' The bivariate exponential geometric distrubution (BEG)
#'
#' Random sample generating function for bivariate exponential geometric distribution with parameters \eqn{\beta > 0} and p in (0,1).
#'
#' rbexpgeo generates random sample from BEG distribution.
#'
#' @param n size of sample.
#' @param beta  numeric paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of random samples generate from BEG model.
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
#' @param beta  numeric paramter which must be numeric greater than  0.
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
  M<-stats:: dgamma(X,shape =N,rate = beta ) *p*(1-p)^N
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
#' @param beta  numeric paramter which must be numeric greater than 0.
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



#' Fits the bivariate exponential geometric distrubution (BEG) to data
#'
#' This function computes the parameter estimates, confidence interval, deviance and covariance matrix for exponential geometric distribution.
#'
#' bexpgeo_fit  fits BEG model to data.
#'
#' @param data is bivariate vector  (X,N) vector representing observations from BEG model
#' @param level confidence level espressed between 0 and 1 (Default is 0.95)
#'
#' @return  list of parameter estimates, confidence interval, deviance and covariance matrix
#'
#' @examples
#' Data.df<-rbexpgeo(n=100,beta = 10,p=0.45)
#' fit<-bexpgeo_fit(Data.df,level = 0.95)
#' fit
#'
#'@references  Kozubowski, T.J., & Panorska, A.K. (2005). A Mixed bivariate distribution with exponential and geometric marginals. Journal of Statistical Planning and Inference, 134, 501-520.
#' \url{https://doi.org/10.1016/j.jspi.2004.04.010}
#'
#' @export
bexpgeo_fit <- function(data,level=0.95) ## data has to be a vector (X,N)
{
  N<-data[,2] ## Get discrete data
  X<-data[,1] ## get continous data
  n<-length(N)
  qt<-(1-level)/2
  qz<- stats:: qnorm(qt)
  z<- abs(qz)
  b<-mean(N)/mean(X)
  p<-1/mean(N)
  sigmabb<-1/(b*b*p)
  sigmapp<-1/((1-p)*p*p)
  J= solve(matrix(c(sigmabb,0,0,sigmapp),byrow = 2, ncol = 2))
  lowerb<-b-z*sqrt(J[1,1]/n)
  lowerp<-p-z*sqrt(J[2,2]/n)
  upperb<-b+z*sqrt(J[1,1]/n)
  upperp<-p+z*sqrt(J[2,2]/n)
  log.like<- sum(log(dbexpgeo(data = data,beta = b,p=p)))
  Output<-data.frame(matrix(c(b,p)),matrix(c(lowerb,lowerp)),matrix(c(upperb,upperp)))
  colnames(Output)<-c("estimate",paste(level*100,"%", " lower bound", sep=""),
                      paste(level*100,"%", " upper bound", sep=""))
  row.names(Output)<- c("beta","p")
  result <- list(Estimates=Output,log.like=log.like,Inverse.Fisher.Matrix=J)
  return(result)
}



