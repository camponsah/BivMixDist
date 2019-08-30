#' The bivariate inverse-Gaussian geometric distribution (BIGG)
#'
#' Random sample generating function for bivariate inverse-Gaussian geometric distribution with parameters \eqn{\mu >0, \phi > 0} and p in (0,1).
#'
#' rbinvgaussgeo generates random samples from BIGG distribution.
#'
#' @param n size of sample.
#' @param mu  shape paramter which must be numeric greater than  0.
#' @param phi  scale paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of samples generate from BGG distribution.
#'
#' @examples
#' data.df<-rbinvgaussgeo(20,mu=1.5, phi=2, p=0.6)
#' data.df
#'
#'@references  Barreto-Souza, W. and Silva R. B. (2019). A bivariate infinitely divisible law for modeling the magnitude and duration of monotone periods of log-returns . Statistica Neerlandica, 73:211-233.
#' \url{https://doi.org/10.1111/stan.12166}
#'
#' @export
rbinvgaussgeo<- function(n,mu,phi,p){
  N<- stats:: rgeom(n,p)+1
  X<- statmod:: rinvgauss(n,mean = N*mu, shape =  phi*N^2)
  return(data.frame(X,N))
}


#' The bivariate inverse-Gaussian geometric distribution (BIGG)
#'
#' Density function for bivariate inverse-Gaussian geometric distribution with parameters \eqn{\mu >0, \phi > 0} and p in (0,1).
#'
#' dbinvgaussgeo is the density function for BGG model.
#'
#' @param data  bivariate vector  (X,N) observations from BGG model.
#' @param mu  shape paramter which must be numeric greater than  0.
#' @param phi  scale paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of densities.
#'
#' @examples
#' data.df<-rbinvgaussgeo(20,mu=1.5, phi=2, p=0.6)
#' den<-rbinvgaussgeo(20,mu=1.5, phi=2, p=0.6)
#' den
#'
#'@references  Barreto-Souza, W. and Silva R. B., (2019). A bivariate infinitely divisible law for modeling the magnitude and duration of monotone periods of log-returns . Statistica Neerlandica, 73:211-233.
#' \url{https://doi.org/10.1111/stan.12166}
#'
#' @export
dbinvgaussgeo<- function(data,mu,phi,p,log.p=FALSE){
  N<-data[,2]
  X<-data[,1]
  M<- statmod:: dinvgauss(X,mean = N*mu, shape =  phi*N^2)*p*((1-p)^(N-1))
  if (log.p == FALSE){
    return(M)
  }else{
    M<- log(M)
    return(M)
  }
}

#' The bivariate inverse-Gaussian geometric distribution (BIGG)
#'
#' Distribution function for bivariate inverse-Gaussian geometric distribution with parameters \eqn{\mu >0, \phi > 0} and p in (0,1).
#'
#' pbgammageo is the distribution function for BGG model.
#'
#' @param data  bivariate vector  (X,N) observations from BGG model.
#' @param mu  shape paramter which must be numeric greater than  0.
#' @param phi  scale paramter which must be numeric greater than  0.
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
#'@references  Barreto-Souza, W. and Silva R. B. (2019). A bivariate infinitely divisible law for modeling the magnitude and duration of monotone periods of log-returns . Statistica Neerlandica, 73:211-233.
#' \url{https://doi.org/10.1111/stan.12166}
#'
#' @export
pbinvgaussgeo<- function(data,mu,phi,p, lower.tail=TRUE,log.p=FALSE){
  cdf<-function(y){
    j <- seq(1:y[2])
    a <- p*((1-p)^j)
    b <-   pnorm(sqrt(phi/y[1])*((y[1]/mu)-j))
    c <- exp(2*phi/mu)*pnorm(-sqrt(phi/y[1])*((y[1]/mu)+j))
    return(a*(b+c))
  }
  M <- apply(data,1, cdf)
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


#' Fits the bivariate inverse-Gaussian geometric distribution (BIGG) to data.
#'
#' This function computes the parameter estimates, confidence interval, deviance and covariance matrix of bivariate inverse-Gaussian geometric distribution.
#''
#' bgammageo_fit BGG model to data.
#'
#' @param data  bivariate vector  (X,N) observations from BGG model.
#' @param level confidence level espressed between 0 and 1 (Default is 0.95).
#'
#' @return  list of parameter estimates, confidence interval, deviance and covariance matrix.
#'
#' @examples
#' Data.df<- rbgammageo(200,alpha=1.5, beta=2, p=0.6)
#' fit <- bgammageo_fit(Data.df)
#' fit
#'
#'@references  Barreto-Souza, W. and Silva R. B. (2019). A bivariate infinitely divisible law for modeling the magnitude and duration of monotone periods of log-returns . Statistica Neerlandica, 73:211-233.
#' \url{https://doi.org/10.1111/stan.12166}
#'
#' @export
bgammageo_fit <- function(data,level=0.95) ## data has to be a vector (X,N)
{
  N<-data[,2]
  X<-data[,1]
  n<-nrow(data)
  qt<-(1-level)/2
  qz<- stats:: qnorm(qt)
  z<- abs(qz)
  mu <- mean(X)/mean(N)
  phi <- n*mu/sum((X-mu*N)^2 /X)
  p <- 1/mean(N)
  kpp <- 1/(p*p*(1-p))
  colnames(J) <- c("mu","phi","p")
  row.names(J) <- c("mu","phi","p")
  lowermu <- mu - z*sqrt((p*mu^3)/(n*phi))
  lowerphi <- phi-z*sqrt(2*phi/n)
  lowerp <- p-z*sqrt(p*p*(1-p)/n)
  uppermu <- mu + z*sqrt((p*mu^3)/(n*phi))
  upperphi <- phi + z*sqrt(2*phi/n)
  upperp <- p + z*sqrt(p*p*(1-p)/n)
  Deviance<- -2*sum(log(dbinvgaussgeo(data,mu,phi,p)))
  Output<-data.frame(matrix(c(mu,phi,p)),matrix(c(lowermu,lowerphi,lowerp)),matrix(c(uppermu,upperphi,upperp)))
  colnames(Output)<-c("estimate",paste(level*100,"%", " lower bound", sep=""),
                      paste(level*100,"%", " upper bound", sep=""))
  row.names(Output)<- c("alpha","beta","p")
  result <- list(Estimates=Output,Deviance=Deviance, Inverse.Fisher.Matrix=J)
  return(result)
}
