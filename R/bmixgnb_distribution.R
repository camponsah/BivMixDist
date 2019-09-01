#' The bivariate mixed gamma negative binomial distribution (BMixGNB)
#'
#' Random sample generating function for bivariate mixed gamma negative binomial distribution with parameters \eqn{\alpha >0, \beta > 0} and p in (0,1).
#'
#' rbgammageo generates random samples from BMixGNB distribution.
#'
#' @param n size of sample.
#' @param alpha numeric paramter which must be numeric greater than  0.
#' @param beta  numeric paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of samples generate from BGG distribution.
#'
#' @examples
#' data.df<-rmixgammanbinom(20,alpha=1.5, beta=2, r=3, p=0.6)
#' data.df
#'
#'@references  Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Lévy process . Journal of Multivariate Analysis, 109:130-145.
#' \url{https://doi.org/10.1016/j.jmva.2012.03.004}
#'
#' @export
rbmixgammanbinom<- function(n,alpha,beta,r,p){
  N<- stats:: rnbinom(n, size = r, p=p)
  X<- stats:: rgamma(n,shape =alpha*(N+r),rate = beta )
  return(data.frame(X,N))
}


#' The bivariate mixed gamma negative binomial distribution  (BMixGNB)
#'
#' Density function for bivariate mixed gamma negative binomial distribution  with parameters \eqn{\alpha >0, \beta > 0} and p in (0,1).
#'
#' dbgammageo is the density function for BMixGNB model.
#'
#' @param data  bivariate vector  (X,N) observations from BMixGNB model.
#' @param alpha numeric paramter which must be numeric greater than  0.
#' @param beta  numeric paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of densities.
#'
#' @examples
#' data.df<-rmixgammanbinom(20,alpha=1.5, beta=2, r=3, p=0.6)
#' den<-dmixgammanbinom(data.df,alpha=1.5, beta=2, r=3, p=0.6)
#' den
#'
#'@references  Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Lévy process . Journal of Multivariate Analysis, 109:130-145.
#' \url{https://doi.org/10.1016/j.jmva.2012.03.004}
#'
#' @export
dbmixgammanbinom<- function(data,alpha,beta,r, p, log.p=FALSE){
  N <- data[,2]
  X <- data[,1]
  M <- (stats:: dgamma(X,shape=alpha*(N+r),rate=beta ))*(stats:: dnbinom(N,size=r,p=p))
  if (log.p == FALSE){
    return(M)
  }else{
    M<- log(M)
    return(M)
  }
}

#' The bivariate mixed gamma negative binomial distribution  (BMixGNB)
#'
#' Distribution function for bivariate mixed gamma negative binomial distribution with parameters \eqn{\alpha >0, \beta > 0} and p in (0,1).
#'
#' pbgammageo is the distribution function for BMixGNB model.
#'
#' @param data  bivariate vector  (X,N) observations from BMixGNB model.
#' @param alpha  numeric paramter which must be numeric greater than  0.
#' @param beta  numeric paramter which must be numeric greater than  0.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x, N \leq n]}, otherwise, \eqn{P[X > x, N > n]}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of distribution.
#'
#' @examples
#' data.df<-rmixgammanbinom(20,alpha=1.5, beta=2, r=3, p=0.6)
#' prob<-pmixgammanbinom(data.df,alpha=1.5, beta=2, r=3, p=0.6)
#' prob
#'
#'@references  Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Lévy process . Journal of Multivariate Analysis, 109:130-145.
#' \url{https://doi.org/10.1016/j.jmva.2012.03.004}
#'
#' @export
pmixgammanbinom<- function(data,alpha,beta,r,p, lower.tail=TRUE,log.p=FALSE){
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


#' Fits the bivariate mixed gamma negative binomial distribution  (BMixGNB) to data.
#'
#' This function computes the parameter estimates, confidence interval, log-likelihood and covariance matrix for the parameters in the bivariate mixed gamma negative binomial distribution
#''
#' bgammageo_fit BMixGNB model to data.
#'
#' @param data  bivariate vector  (X,N) observations from BMixGNB model.
#' @param level confidence level espressed between 0 and 1 (Default is 0.95).
#'
#' @return  list of parameter estimates, confidence interval, deviance and covariance matrix.
#'
#' @examples
#' Data.df<- rbgammageo(200,alpha=1.5, beta=2, p=0.6)
#' fit <- bgammageo_fit(Data.df)
#' fit
#'
#'@references  Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Lévy process . Journal of Multivariate Analysis, 109:130-145.
#' \url{https://doi.org/10.1016/j.jmva.2012.03.004}
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
  ### Log likelihood function
  log.lik <- function(par) { #par[1]=alpha, par[2]=beta
    #b=par*mean(N)/mean(X)
    B<-par*mean(N)/mean(X)
    ll<- par*mean(N)*log(B+1e-16)-B*mean(X)+mean(par*N*log(X)-lgamma(par*N))
    return(-ll)
  }
  ## Estimate alpha
  V<- stats:: var(X)
  alpha<- mean(N*((mean(X))^2)/V)
  constant<- log(mean(N)/mean(X)) + mean(N*log(X/N))/mean(N)
  if (constant<0){
    a<- stats:: nlm(log.lik,p=alpha)$estimate
  }else{
    a<- Inf
    message("MLE of alpha does not exist!")
  }
  b<-a*mean(N)/mean(X) ## estimate beta
  p<-1/mean(N) # estimate p
  ### Sum to infinity function
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
  log.like<-N*a*log(b)-lgamma(a*N)+(a*N-1)*log(X)-b*X+log(p)+(N-1)*log(1-p)
  Deviance<- -2*sum(log.like)
  Output<-data.frame(matrix(c(a,b,p)),matrix(c(lowera,lowerb,lowerp)),matrix(c(uppera,upperb,upperp)))
  colnames(Output)<-c("estimate",paste(level*100,"%", " lower bound", sep=""),
                      paste(level*100,"%", " upper bound", sep=""))
  row.names(Output)<- c("alpha","beta","p")
  result <- list(Estimates=Output,Deviance=Deviance, Inverse.Fisher.Matrix=J)
  return(result)
}
