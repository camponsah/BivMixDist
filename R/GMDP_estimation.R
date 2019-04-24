#' Fits the gamma mixture distcrete Pareto distribution (GMDP) to data.
#'
#' This function computes the estimates  of parameter in GMDP model.
#'
#' bgammageo_fit GMDP model to data.
#'
#' @param data  bivariate vector  (X,N) observations from GMDP model.
#' @param delta  initial guess of true parameter \eqn{\delta} which is numeric and must be greater than  0 (Default value is 1).
#' @param p  nitial guess of true parameter p which must be numeric value between 0 and 1 (Default value is 0.5).
#' @param method method of estimation: eqn{EM=}EM algorithm or eqn{MLE=}maximum likelihood estimation (Default method is EM).
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
gammamixdpareto_fit <- function(data,delta=1, p=0.5, method="EM") ## data has to be a vector (X,N)
{
  N<-data[,2]
  X<-data[,1]
  n<-nrow(data)
  #qt<-(1-level)/2
  #qz<- stats:: qnorm(qt)
  #z<- abs(qz)
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
    message("MLE of alpha does not exist!")
  }
  b <- a*mean(N)/mean(X) ## estimate beta
  if(method == "EM"){
    fit <- dpareto_em(N, maxiter = 1000)
    delta <- fit$par$delta
    p <- fit$par$p
  }else{
    fit <- stats:: optim(par = c(delta, p), fn = log.lik.DP)
    delta <- fit$par[1]
    p <- fit$par[2]
  }

  Output<-data.frame(t(matrix(c(a,b,delta,p))))
  colnames(Output)<-c("alpha"," beta", "delta", "p")
  return(Output)
}
