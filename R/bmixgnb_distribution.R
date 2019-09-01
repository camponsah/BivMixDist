#' The bivariate mixed gamma negative binomial distribution (BMixGNB)
#'
#' Random sample generating function for bivariate mixed gamma negative binomial distribution with parameters \eqn{\alpha >0, \beta > 0} and p in (0,1).
#'
#' rbmixgammanbinom generates random samples from BMixGNB distribution.
#'
#' @param n size of sample.
#' @param alpha numeric paramter which must be numeric greater than  0.
#' @param beta  numeric paramter which must be numeric greater than  0.
#' @param r  integer paramter which must be greater than or equal to  0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of samples generate from BGG distribution.
#'
#' @examples
#' data.df<-rbmixgammanbinom(20,alpha=1.5, beta=2, r=3, p=0.6)
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
#' dbmixgammanbinom is the density function for BMixGNB model.
#'
#' @param data  bivariate vector  (X,N) observations from BMixGNB model.
#' @param alpha numeric paramter which must be numeric greater than  0.
#' @param beta  numeric paramter which must be numeric greater than  0.
#' @param r  integer paramter which must be greater than or equal to  0.
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
#' pmixgammanbinom is the distribution function for BMixGNB model.
#'
#' @param data  bivariate vector  (X,N) observations from BMixGNB model.
#' @param alpha  numeric paramter which must be numeric greater than  0.
#' @param beta  numeric paramter which must be numeric greater than  0.
#' @param r  integer paramter which must be greater than or equal to  0.
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
pmixgammanbinom <- function(data,alpha,beta,r,p, lower.tail=TRUE,log.p=FALSE){
  cdf <- function(y){
    j <- seq(0:y[2])
    S <-  ((1-p)^j)*(gamma(j+r)/gamma(j+1)) * pracma::gammainc(beta*y[1],alpha*(r+j))[3]
    return((p^r/gamma(r))*sum(S))
  }
  M <- apply(data, 1 , cdf)
  if (lower.tail==FALSE & log.p==TRUE){
    M <- log(1-M)
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
#' bmixgammanbinom_fit BMixGNB model to data.
#'
#' @param data  bivariate vector  (X,N) observations from BMixGNB model.
#' @param r  numeric parameter, which must NULL if the parameter in unknown or an integer.
#' @param level confidence level espressed between 0 and 1 (Default is 0.95).
#'
#' @return  list of parameter estimates, confidence interval, deviance and covariance matrix.
#'
#' @examples
#' data.df<-rbmixgammanbinom(20,alpha=1.5, beta=2, r=3, p=0.6)
#' fit <- bmixgammanbinom_fit(data.df)
#' fit
#'
#'@references  Barreto-Souza, W. (2012). Bivariate gamma-geometric law and its induced Lévy process . Journal of Multivariate Analysis, 109:130-145.
#' \url{https://doi.org/10.1016/j.jmva.2012.03.004}
#'
#' @export
bmixgammanbinom_fit <- function(data,r=NULL,level=0.95) ## data has to be a vector (X,N)
{
  N <- data[,2]
  X <- data[,1]
  n <- length(N)
  qt <- (1-level)/2
  qz <- stats:: qnorm(qt)
  z <- abs(qz)
  log.lik1 <- function(par) {
    alpha <- par[1]
    beta <- par[2]
    r <- par[3]
    p <- par[4]
    ll <- dbmixgammanbinom(data,alpha,beta,r, p, log.p=TRUE)
    return(-sum(ll))
  }
  log.lik2 <- function(par) {
    alpha <- par[1]
    beta <- par[2]
    p <- par[3]
    ll <- dbmixgammanbinom(data,alpha,r, p, log.p=TRUE)
      return(-sum(ll))
  }
  sumToInfinityaa <- function(alpha,beta,r,p){
    j <- 0
    error <- 1
    S1 <- 0
    while(error > 0.00001){
      S <- S1+ ((r+j)^2)*(1-p)^j* psigamma(alpha*(j+r), deriv = 1) *gamma(r+j)/gamma(j+1)
      j <- j+1
      error <- abs(S-S1)
      S1 <- S
    }
    return((p^r /gamma(r))*S1)
  }
  sumToInfinityar <- function(alpha,r,p){
    j <- 0
    error <- 1
    S1 <- 0
    while(error > 0.00001){
      S <- S1+ (r+j)*(1-p)^j* psigamma(alpha*(j+r), deriv = 1) *gamma(r+j)/gamma(j+1)
      j <- j+1
      error <- abs(S-S1)
      S1 <- S
    }
    return((alpha*p^r /gamma(r))*S1)
  }
  sumToInfinityrr <- function(alpha,r,p){
    j <- 0
    error <- 1
    S1 <- 0
    while(error > 0.00000001){
      t <- (1-p)^j * gamma(r+j)/gamma(j+1)
      S <- S1 + t* (alpha^2 *psigamma(alpha*(j+r), deriv = 1) -psigamma((j+r), deriv = 1) )
      j <- j+1
      error <- abs(S-S1)
      S1 <- S
    }
    return(psigamma(r, deriv = 1) + (p^r /gamma(r))*S1)
  }
  if(!is.null(r)){
    par <- c(alpha, beta,p)
    par<- try(suppressWarnings(nlm(f=log.lik2,p=par, ndigit = 12)$estimate),
              silent=TRUE)
    alpha <- par[1]
    beta <- par[2]
    p <- par[3]
    kaa <- sumToInfinityaa(alpha = alpha, r=r, p=p)
    kbb <- (alpha*r)/(p*beta^2)
    kpp <- r/(p*p(1-p))
    kba <- -r/(p*beta)
    J <- rbind(c(kaa,kba,0), c(kba,kbb,0), c(0,0,kpp))
    J <- solve(J)
    J <- data.frame(J)
    colnames(J) <- c("alpha","beta","p")
    row.names(J)<- c("alpha","beta","p")
    lowa <- alpha-z*sqrt(J[1,1]/n)
    lowb <- beta-z*sqrt(J[2,2]/n)
    lowp <- p - z*sqrt(J[3,3]/n)
    upa <- alpha + z*sqrt(J[1,1]/n)
    upb <- beta + z*sqrt(J[2,2]/n)
    upp <- p + z*sqrt(J[3,3]/n)
    log.like <- -log.lik2(par = c(alpha,beta,r,p))
    Output <- rbind(c(alpha,beta,p),c(lowa,lowb,lowp),c(upa,upb, upp))
colnames(Output) <- c("estimate",paste(level*100,"%", " lower bound", sep=""),
                      paste(level*100,"%", " upper bound", sep=""))
row.names(Output) <- c("alpha","beta","p")
  }
  else{
    par <- c(alpha, beta, r,p)
    par <- try(suppressWarnings(nlm(f=log.lik1,p=par, ndigit = 12)$estimate),
               silent=TRUE)
    alpha <- par[1]
    beta <- par[2]
    r < - par[3]
    p <- par[4]
    kaa <- sumToInfinityaa(alpha = alpha, r=r, p=p)
    kbb <- (alpha*r)/(p*beta^2)
    krr <- sumToInfinityrr(alpha = alpha, r=r, p=p)
    kpp <- r/(p*p(1-p))
    kba <- -r/(p*beta)
    kar <- sumToInfinityar(alpha = alpha, r=r, p=p)
    kbr <- -alpha/beta
    kpr <- -1/p
    J <- rbind(c(kaa,kba,kar,0), c(kba,kbb,kbr,0),
               c(kar,kbr,krr,kpr), c(0,0,kpr,kpp))
    J<- solve(J)
    J<-data.frame(J)
    colnames(J) <- c("alpha","beta","r","p")
    row.names(J) <- c("alpha","beta","r","p")
    lowa <- alpha - z*sqrt(J[1,1]/n)
    lowb <- beta - z*sqrt(J[2,2]/n)
    lowr <- r - z*sqrt(J[3,3]/n)
    lowp <- p - z*sqrt(J[4,4]/n)
    upa <- alpha + z*sqrt(J[1,1]/n)
    upb <- beta + z*sqrt(J[2,2]/n)
    upr <- r + z*sqrt(J[3,3]/n)
    upp <- p - z*sqrt(J[4,4]/n)
    log.like <- -log.lik1(par = c(alpha,beta,r,p))
    Output <- rbind(c(alpha,beta,r,p),c(lowa,lowb,lowr,lowp),c(upa,upb,upr, upp))
    colnames(Output) <- c("estimate",paste(level*100,"%", " lower bound", sep=""),
                          paste(level*100,"%", " upper bound", sep=""))
    row.names(Output) <- c("alpha","beta","r","p")
  }
  result <- list(Estimates=Output,log.like=log.like, Inverse.Fisher.Matrix=J)
  return(result)
}
