#' The gamma mixture discrete Pareto distribution (GMDP)
#'
#' Random sample generating function for gamma mixture discrete Pareto distribution with parameters \eqn{\alpha >0, \beta > 0, \delta \ge 0} and p in (0,1).
#'
#' rgammamixdpareto generates random samples from GMDP distribution.
#'
#' @param n size of sample.
#' @param alpha  numeric parameter which must be greater than  0.
#' @param beta  numeric parameter which must be greater than  0.
#' @param  delta numeric parameter which must be greater than or equal to  0.
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
  N <- rdpareto(n, delta = delta,p=p)
  X<- stats:: rgamma(n,shape =alpha*N,rate = beta )
  return(data.frame(X,N))
}

#' The gamma mixture discrete Pareto distribution (GMDP)
#'
#' Density function for gamma mixture discrete Pareto distribution with parameters \eqn{\alpha >0, \beta > 0, \delta \ge 0} and p in (0,1).
#'
#' dgammamixdpareto is the density function of GMDP model.
#'
#' @param data  vector  (X,N) representing observation from GMDP model.
#' @param alpha  numeric parameter which must be greater than  0.
#' @param beta  numeric parameter which must be greater than  0.
#' @param  delta numeric parameter which must be greater than or equal to  0.
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
#'@references  Amponsah, C. K., Kozubowski, T. J. and Panorska, A. K. (2019). A bivariate distribution with gamma mixture discrete Pareto margins . In print.
#' \url{https://scholarworks.unr.edu/bitstream/handle/11714/2065/Amponsah_unr_0139M_12378.pdf?sequence=1&isAllowed=y}
#'
#' @export
dgammamixdpareto<- function(data,alpha,beta,delta,p,log.p=FALSE){
  if (is.numeric(data[,1]) && is.numeric(data[,2])){
    N <- as.integer(data[,2])
    X<- as.double(data[,1])
  }
  else stop("all entries of argument 'data' must be numeric")
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

#' The gamma mixture discrete Pareto distribution (GMDP)
#'
#' Distribution function for gamma mixture discrete Pareto distribution with parameters \eqn{\alpha >0, \beta > 0, \delta \ge 0} and p in (0,1).
#'
#' pgammamixdpareto is the distribution function of GMDP model.
#'
#' @param data  vector  (X,N) representing observation from GMDP model.
#' @param alpha  numeric parameter which must be greater than  0.
#' @param beta  numeric parameter which must be greater than  0.
#' @param  delta numeric parameter which must be greater than or equal to  0.
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
  if (is.numeric(data[,1]) && is.numeric(data[,2])){
    N <- as.integer(data[,2])
    X<- as.double(data[,1])
  }
  else stop("all entries of argument 'data' must be numeric")
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


#' Fits the gamma mixture discrete Pareto distribution (GMDP) to data.
#'
#' This function computes the estimates  of parameter in GMDP model.
#'
#' bgammageo_fit GMDP model to data.
#'
#' @param data  bivariate vector  (X,N) observations from GMDP model.
#' @param delta  initial guess of true parameter \eqn{\delta} which is numeric and must be greater than  0 (Default value is 1).
#' @param p  initial guess of true parameter p which must be numeric value between 0 and 1 (Default value is 0.5).
#' @param method method of estimation for discrete Pareto parameter: \eqn{EM=}EM algorithm or \eqn{MLE=}maximum likelihood estimation (Default method is EM).
#'
#' @return  vector of parameter estimates.
#'
#' @examples
#' Data.df<- rgammamixdpareto(n=1000,alpha=1,beta=2,delta=0.45, p=0.8)
#' fit <- gammamixdpareto_fit(data = Data.df)
#' fit
#'
#'@references  Amponsah, C. K., Kozubowski, T. J. and Panorska, A. K. (2019). A bivariate distribution with gamma mixture discrete Pareto marginals . In print.
#' \url{https://scholarworks.unr.edu/bitstream/handle/11714/2065/Amponsah_unr_0139M_12378.pdf?sequence=1&isAllowed=y}
#'
#' @export
gammamixdpareto_fit <- function(data,delta=1, p=0.5, method="MLE")
{
  if (is.numeric(data[,1]) && is.numeric(data[,2])){
    N <- as.integer(data[,2])
    X<- as.double(data[,1])
  }
  else stop("all entries of argument 'data' must be numeric")
  method <- match.arg(method)
  n <- length(N)
  ### Log likelihood function
  log.lik <- function(par) {
    B<-par*mean(N)/mean(X)
    ll<- par*mean(N)*log(B+1e-16)-B*mean(X)+mean(par*N*log(X)-lgamma(par*N))
    return(-ll)
  }
  log.lik.DP <- function(par) { #par[2]=p and par[1]=delta
    ll1<- (1/(1-par[1]*(N-1)*log(1-par[2])))^(1/par[1])
    ll2<- (1/(1-par[1]*N*log(1-par[2])))^(1/par[1])
    D<- sum(log(ll1-ll2))
    -D
  }
  ## Estimate alpha
  V <- stats::var(X)
  alpha<- mean(N*((mean(X))^2)/V)
  constant<- log(mean(N)/mean(X)) + mean(N*log(X/N))/mean(N)
  if (constant<0){
    a<- stats:: nlm(log.lik,p=alpha)$estimate
  }
  else{
    a<- Inf
    #message("MLE of alpha does not exist!")
  }
  b <- a*mean(N)/mean(X) ## estimate beta
  if(method == "MLE"){
    fit <- stats:: nlm(log.lik.DP,p=c(delta, p))
    delta <- fit$estimate[1]
    p <- fit$estimate[2]
  }else{
    fit <- dpareto_em(N)
    delta <- fit$par$delta
    p <- fit$par$p
  }

  Output<-data.frame(t(matrix(c(a,b,delta,p))))
  colnames(Output)<-c("alpha"," beta", "delta", "p")
  return(Output)
}


