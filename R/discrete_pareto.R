#' The discrete Pareto distribution
#'
#' Random sample generating function for discrete Pareto distribution with parameters \eqn{\delta \ge 0} and p in (0,1).
#'
#' rdpareto generates random sample from discrete Pareto distribution.
#'
#' @param n size of sample.
#' @param delta  shape parameter which must be numeric greater than or equal to 0.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of samples generate from discrete Pareto distribution.
#'
#' @examples
#' N<-rdpareto(20, delta=0.2, p=0.3)
#' N
#'
#'@references  Buddana, A. and  Kozubowski, T. J. (2014). Discrete Pareto distribution. Journal of Economics and Quality Control, 29(2):143-156.
#' \url{https://doi.org/10.1515/eqc-2014-0014}
#'
#' @export
rdpareto<-function(n,delta,p){
  u<- stats:: runif(n)
  sigma<- - 1/(delta*log(1 - p))
  N <- ceiling(sigma*((1 - u)^(- delta) - 1))
  return(N)
}

#' The discrete Pareto distribution
#'
#' Probability mass function for discrete Pareto distribution with parameters \eqn{\delta \ge 0} and p in (0,1).
#'
#' ddpareto gives the probability mass.
#'
#' @param N A vector of random sample from discrete Pareto.
#' @param delta A shape parameter which must be numeric greater than or equal to 0.
#' @param p A numeric parameter between 0 and 1
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return vector representing probabilities \eqn{P(N=n)}.
#'
#' @examples
#' prob<-ddpareto(seq(1,10,1),delta=0.2,p=0.6)
#' prob
#'
#'@references  Buddana, A. and  Kozubowski, T. J. (2014). Discrete Pareto distribution. Journal of Economics and Quality Control, 29(2):143-156.
#' \url{https://doi.org/10.1515/eqc-2014-0014}
#'
#' @export
ddpareto<-function(N,delta,p,log.p=FALSE){
 W1<- (1 - delta*(N - 1)*log(1 - p))^( - 1/delta)
 W2<- (1 - delta*N*log(1 - p))^( - 1/delta)
 M<- W1 - W2
 if (log.p == FALSE){
   return(M)
   }else{
   M<- log(M)
   return(M)
   }
}

#' The discrete Pareto distribution
#'
#' Distribution function for discrete Pareto distribution with parameters \eqn{\delta \ge 0} and p in (0,1).
#'
#' pdpareto gives the distribution function.
#'
#' @param q vector of quantiles representing numbers from discrete Pareto distribution.
#' @param delta  shape parameter which must be numeric greater than or equal to 0.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[N\leq n]}, otherwise, \eqn{$P[N> n]$}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return A vector cumulative probabilities
#'
#' @examples
#' prob<-pdpareto(seq(1,10, 1),delta=0.2,p=0.6)
#' prob
#'
#'@references  Buddana, A. and  Kozubowski, T. J. (2014). Discrete Pareto distribution. Journal of Economics and Quality Control, 29(2):143-156.
#' \url{https://doi.org/10.1515/eqc-2014-0014}
#'
#' @export
pdpareto<- function(q,delta,p,lower.tail=TRUE,log.p=FALSE){
  sigma = -1/(delta*log(1 - p))
  M<- 1- ((1 + q/sigma)^(- 1/delta))
if (lower.tail == TRUE & log.p == FALSE){
  return(M)
  }else if (lower.tail == TRUE & log.p == TRUE){
  M<- log(M)
  return(M)
  }else if (lower.tail == FALSE & log.p == TRUE){
  M<- log(1 - M)
  return(M)
  }else {
  return(1 - M)
  }
}

#' The discrete Pareto distribution
#'
#' Quantile function for discrete Pareto distribution with parameters \eqn{\delta \ge 0} and p in (0,1).
#'
#' qdpareto gives the quantile function.
#'
#' @param prob Vector of probabilities.
#' @param delta  shape parameter which must be numeric greater than or equal to 0.
#' @param p  numeric parameter between 0 and 1.
#'
#' @return  vector quantiles from discrete Pareto distribution.
#'
#' @examples
#' q<-rdpareto(seq(0.1,0.6,0.1), delta=0.2,p=0.6)
#' q
#'
#'@references  Buddana, A. and  Kozubowski, T. J. (2014). Discrete Pareto distribution. Journal of Economics and Quality Control, 29(2):143-156.
#' \url{https://doi.org/10.1515/eqc-2014-0014}
#'
#' @export
qdpareto<-function(prob,delta,p){
  sigma <- - 1/(delta*log(1 - p))
  return(ceiling(sigma*((1-  prob)^(- delta) - 1)))
}



#' EM algorithm function for discrete Pareto distribution
#'
#' This function computes the parameter estimates of discrete Pareto distribution using the EM algorithm.
#'
#' Takes initial guess for the parameters in discrete Pareto and the algorithm will estimate the MLE.
#'
#' @param N  vector of random sample from discrete Pareto distribution
#' @param delta shape parameter which must be numeric greater than or equal to 0
#' @param p  numeric parameter between 0 and 1
#' @param maxiter maximum number of iterations
#' @param tol tolerance value
#' @param verb If TRUE, estimates are printed during each iteration.
#'
#' @return  list containg parameter estimates, Deviance and data frame of iteration
#'
#' @examples
#' N<-rdpareto(500, delta=0.2,p=0.6)
#' fit<-dpareto_em(N,maxiter = 1000)
#' fit$par
#'
#'@references  Amponsah, C. K.,  Kozubowski, T. J. and Panorska (2019). A computational approach to estimation of discrete Pareto parameters. Inprint.
#'
#' @export
dpareto_em <- function(N, delta = 1, p = NULL, maxiter = 1000,
                       tol = 1e-8, verb=FALSE){
  n <- length(N)
  if (!is.numeric(delta)) stop("argument 'delta' must be numeric greater than 0")
  if (is.null(p)){
    p <- 1/mean(N)
  }
  if ((p >= 1) | (p <= 0)) stop("argument 'p' must be numeric between 0 and 1")
  gamma_1<- - 1/(delta*log(1 - p))
  eta<- 1/delta
  func_eta<-function(eta){
    ll<- eta*log(eta) - lgamma(eta)- eta - eta*log(mean(a)) + eta*mean(c)
    return( -ll)
  }
  ll_old <- sum(log(ddpareto(N=N,delta=delta, p=p)))
  k = 0
  output <- c(k,delta, p, ll_old)
  diff <- tol +1
  par <- c(delta,p)
  while(diff > tol & k < maxiter){
    ### E step
    const<- 1/(gamma(eta)*((gamma_1 + N - 1)^(-eta) - (gamma_1 + N)^(-eta)))
    a<- const * gamma(eta + 1) * ((gamma_1 + N - 1)^(-(eta + 1)) - (gamma_1 + N)^(-(eta + 1)))
    t1<- digamma(eta) - log(gamma_1 + N - 1)
    t2<- digamma(eta) - log(gamma_1 + N)
    c<-  const*gamma(eta) * (t1*(gamma_1 + N - 1)^(-eta) - t2*(gamma_1 + N)^(-eta))
    #### M step
    constant<-mean(c)-log(mean(a))
    if(constant<0){
      eta=try(suppressWarnings(stats:: nlm(f=func_eta,p=eta, ndigit = 12)$estimate),
              silent=TRUE)
    }
    else{
      eta <- Inf
    }
    delta<- 1/eta
    p <- 1 - exp( - mean(a))
    gamma_1 <- - 1/(delta*log(1 - p))
    ll_new <- sum(log(ddpareto(N=N,delta=delta, p=p)))
    diff <- (sum((c(delta,p)-par)^2))^(1/2)
    par <- c(delta,p)
    ll_old <- ll_new
    k <- k + 1
    output <- rbind(output, c(k, delta, p, ll_new))
    if (verb) {
      cat("iteration =", k, " log.lik.diff =", diff, " log.lik =",
          ll_new, "\n")
    }
  }
  if (k == maxiter) {
    cat("Warning! Convergence not achieved!", "\n")
  }
  output <- data.frame(output)
  colnames(output) <- c("iteration","delta","p","log-lik")
  par<-data.frame(t(c(delta,p)))
  colnames(par)<-c("delta","p")
  result <- list(par=par, log.like=ll_new, iterations=k, output=output)
  return(result)
}


#' Fisher Scoring algorithm (FSA) for estimation of discrete Pareto distribution parameters.
#'
#' This function computes MLE of discrete Pareto distribution using the FSA algorithm.
#'
#' Takes initial guess for the parameters in discrete Pareto and the algorithm will estimate the MLE.
#'
#' @param N  vector of random sample from discrete Pareto distribution
#' @param delta shape parameter which must be numeric greater than or equal to 0
#' @param p  numeric parameter between 0 and 1
#' @param lambda  learning parameter between 0 and 1
#' @param maxiter maximum number of iterations
#' @param tol tolerance value
#' @param verb If TRUE, estimates are printed during each iteration.
#'
#' @return  list containg parameter estimates, Deviance and data frame of iteration
#'
#' @examples
#' N<-rdpareto(500, delta=0.2,p=0.6)
#' fit<-dpareto_em(N,maxiter = 1000)
#' fit$par
#'
#'@references  Amponsah, C. K.,  Kozubowski, T. J. and Panorska (2019). A computational approach to estimation of discrete Pareto parameters. Inprint.
#'
#' @export
dpareto_fsa <- function(N, delta = 1, p = NULL, lambda=NULL, maxiter = 20,
                       tol = 1e-8, verb=FALSE){
  #n <- length(N)
  if (!is.numeric(delta)) stop("argument 'delta' must be numeric greater than 0")
  if (is.null(p)){
    p <- 1/mean(N)
  }
  if ((p >= 1) | (p <= 0)) stop("argument 'p' must be numeric between 0 and 1")
  h <- function(delta, p,n){
    a <- log(1-delta*n*log(1-p))/(delta^2)
    b <- (n*log(1-p))/(delta*(1-delta*n*log(1-p)))
    c <- (1-delta*n*log(1-p))^(-1/delta)
    return((a+b)*c)
  }
  w <- function(delta, p,n){
    a <- (n/(1-p))*(1-delta*n*log(1-p))^(-(1+ 1/delta))
   return(a)
  }
  g <- function(delta, p,n){
    a <- log(1-delta*n*log(1-p))
    b <- (1-delta*n*log(1-p))^(1/delta)
    return(a/b)
  }
  v <- function(delta, p,n){
    b <- (1-delta*n*log(1-p))^(1 + 1/delta)
    return(n/b)
  }
  score <- function(delta,p,N){
    dL <- (h(delta,p,N-1)-h(delta,p,N))/ddpareto(N,delta,p)
    pL <- (w(delta,p,N-1)-w(delta,p,N))/ddpareto(N,delta,p)
    return(matrix(c(sum(dL),sum(pL)),nrow = 1, ncol = 2))
  }
  I_mat <- function(delta, p, N){
    D <- matrix(c(-1/(delta^2), -log(1-p), 0, delta/(1-p)),nrow = 2, ncol = 2)
    a <- g(delta,p,N-1)-g(delta,p,N)
    b <- v(delta,p,N-1)-v(delta,p,N)
    f <- ddpareto(N,delta,p)
    w11 <- mean((a/f)^2)
    w22 <- mean((b/f)^2)/(delta^2)
    w12 <- mean((a*b/(f^2)))/delta
    A <- matrix(c(w11,w12,w12,w22),nrow = 2, ncol = 2)
    return((t(D)%*%A)%*%D)
  }
  log_like<- function(delta, p){
    W1<- (1 - delta*(N - 1)*log(1 - p))^( - 1/delta)
    W2<- (1 - delta*N*log(1 - p))^( - 1/delta)
    return(sum(log( W1 - W2)))
  }
  ll_old <- log_like(delta, p)
   theta <-  c(delta,p)
  k = 0
  output <- c(k,delta, p, ll_old)
  diff <- tol +1
  par <- c(delta,p)
  while(diff > tol && k < maxiter ){ #diff > tol && & k < maxiter
     S <- matrix(c(score(par[1],par[2],N)), nrow = 2,ncol = 1)
     I <- I_mat(par[1],par[2],N)
    par <- par + lambda* t(S) %*% solve(I+(S%*%t(S))*diag(2)) # +0.0001*diag(2)
    ll_new <- log_like(par[1], par[2])
    diff <- sqrt(sum((par - theta)^2))
    theta <- par
    k <- k + 1
    output <- rbind(output, c(k, par[1], par[2], ll_new))
    if (verb) {
      cat("iteration =", k, " log.lik.diff =", diff, " log.lik =",
          ll_new, "\n")
    }
  }
  if (k == maxiter) {
    cat("Warning! Convergence not achieved!", "\n")
  }
  output <- data.frame(output)
  colnames(output) <- c("iteration","delta","p","log-lik")
  par<-data.frame(t(c(par[1],par[2])))
  colnames(par)<-c("delta","p")
  result <- list(par=par, log.like=ll_new, iterations=k, output=output)
  return(result)
}

