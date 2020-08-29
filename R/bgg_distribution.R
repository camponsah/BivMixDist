#' The bivariate gamma geometric distribution (BGG)
#'
#' Random sample generating function for bivariate gamma geometric distribution with parameters \eqn{\alpha >0, \beta > 0} and p in (0,1).
#'
#' rbgammageo generates random samples from BGG distribution.
#'
#' @param n size of sample.
#' @param alpha  numeric parameter which must be greater than  0.
#' @param beta  numeric parameter which must be  greater than  0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of random samples generate from BGG distribution.
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
#' @param alpha  numeric parameter which must be  greater than  0.
#' @param beta  numeric parameter which must be  greater than  0.
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
  N <- data[,2]
  X <- data[,1]
  M <- ((1-p)^(N-1))*(p*beta^(alpha*N))*(X^(alpha*N-1))*exp(-beta*X)/gamma(alpha*N)
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
#' @param alpha  numeric parameter which must be  greater than  0.
#' @param beta  numeric parameter which must be  greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x, N \leq n]}, otherwise, \eqn{P[X > x, N > n]}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of distribution.
#'
#' @examples
#' data.df <- rbgammageo(20, alpha=1.5, beta=2, p=0.6)
#' prob <- pbgammageo(data.df, alpha=1.5, beta=2, p=0.6)
#' prob
#'
#'@references  Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Lévy process . Journal of Multivariate Analysis, 109:130-145.
#' \url{https://doi.org/10.1016/j.jmva.2012.03.004}
#'
#' @export
pbgammageo<- function(data, alpha, beta, p, lower.tail=TRUE, log.p=FALSE){
  cdf <- function(y) {
    j <- seq(1,y[2])
    for (i in 1:length(j)){
      S <- pracma:: gammainc(beta*y[1],j[i]*alpha)[3] * p*((1-p)^(j-1))
      return(sum(S))
    }
  }
  M <- apply(data,1,cdf)
  if (lower.tail==FALSE & log.p==TRUE){
    return(log(1-M))
  }
  else if (lower.tail==TRUE & log.p==TRUE){
    return(log(M))
  }
  else if (lower.tail==TRUE & log.p==FALSE){
    return(M)
  }
  else {
    return(1-M)
  }
}


#' Fits the bivariate gamma geometric distribution (BGG) to data.
#'
#' This function computes the parameter estimates, confidence interval, deviance and covariance matrix of gamma geometric distribution.
#''
#' bgammageo_fit BGG model to data.
#'
#' @param data  bivariate vector  (X,N) observations from BGG model.
#' @param level confidence level expressed between 0 and 1 (Default is 0.95).
#'
#' @return  list of parameter estimates, confidence interval, deviance and covariance matrix.
#'
#' @examples
#' Data.df<- rbgammageo(200,alpha=1.5, beta=2, p=0.6)
#' fit <- bgammageo_fit(Data.df)
#' fit
#'
#'@references  Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Levy process . Journal of Multivariate Analysis, 109:130-145.
#' \url{https://doi.org/10.1016/j.jmva.2012.03.004}
#'
#' @export
bgammageo_fit <- function(data,level=0.95) ## data has to be a vector (X,N)
{
  N<-data[,2]
  X<-data[,1]
  n<-length(N)
  qt<-(1-level)/2
  qz<- stats:: qnorm(qt)
  z<- abs(qz)
  log.lik <- function(par) { #par[1]=alpha, par[2]=beta
    B<-par*mean(N)/mean(X)
    ll<- par*mean(N)*log(B+1e-16)-B*mean(X)+mean(par*N*log(X)-lgamma(par*N))
    return(-ll)
  }
  sumToInfinity<- function(p,a){
    j=1
    error=1
    S1=0
    while(error>0.00001){
      S=S1+ p*j*j *((1-p)^(j-1))* psigamma(j*a, deriv = 1)
      j=j+1
      error=abs(S-S1)
      S1=S
    }
    return(S1)
  }
  V<- stats:: var(X)
  alpha<- mean(N*((mean(X))^2)/V)
  constant<- log(mean(N)/mean(X)) + mean(N*log(X/N))/mean(N)
  if (constant<0){
    a <- stats:: nlm(log.lik,p=alpha)$estimate
    b <- a*(mean(N)/mean(X))
    p <- 1/mean(N)
    sigmaaa<- sumToInfinity(p, a)
    sigmabb<- a/(p*b*b)
    sigmaab<- -1/(p*b)
    sigmapp<-1/(p*p*(1-p))
    J= solve(matrix(c(sigmaaa,sigmaab,0,sigmaab,sigmabb,0,0,0,sigmapp),byrow = 3, ncol = 3))
    J<-data.frame(J)
    colnames(J)<- c("alpha","beta","p")
    row.names(J)<- c("alpha","beta","p")
    lowera<-a-z*sqrt(J[1,1]/n)
    lowerb<-b-z*sqrt(J[2,2]/n)
    lowerp<-p-z*sqrt(J[3,3]/n)
    uppera<-a+z*sqrt(J[1,1]/n)
    upperb<-b+z*sqrt(J[2,2]/n)
    upperp<-p+z*sqrt(J[3,3]/n)
    log.like <- sum(log(dbgammageo(data = data,alpha = a, beta = b, p=p)))
    Output<-rbind(c(a,lowera,uppera),c(b,lowerb,upperb),c(p,lowerp,upperp))
    colnames(Output)<-c("estimate",paste(level*100,"%", " lower bound", sep=""),
                        paste(level*100,"%", " upper bound", sep=""))
    row.names(Output)<- c("alpha","beta","p")
    result <- list(estimates=Output,log.like=log.like, Inverse.Fisher.Matrix=J)
    return(result)
  }else{
    cat("Warning! MLE of alpha and beta do not exist!", "\n")
  }

}
