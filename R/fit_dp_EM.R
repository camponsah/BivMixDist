#' EM algorithm function for discrete Pareto distribution
#'
#' This function computes the parameter estimates for exponential geometric distribution using the EM algorithm.
#'
#' Takes initial guess of parameters in discrete pareto and use Em algorithm to estmate the MLE.
#'
#' @param data  vector of random sample from discrete pareto distribution
#' @param delta shape paramter which must be numeric greater than or equal to 0
#' @param p  numeric parameter between 0 and 1
#' @param maxiter maximum number of iterations
#' @param tol tolerance value
#'
#' @return  list containg parameter estimates, Deviance and data frame of iteration
#'
#' @examples
#' N<-rdpareto(500, delta=0.2,p=0.6)
#' fit<-dpareto_em(N,maxiter = 1000)
#' fit$par
#'
#'@references  Amponsah, C. K.,  Kozubowski, T. J. and Panorska (2019). Fitting Discrete Pareto Distribution to data using the EM Algorithm. Unpublished.
#'
#' @export
dpareto_em <- function(data, delta = 1, p = 0.5, maxiter = 500, tol = 1e-16){
  N<- data
  n <- length(N)
  gamma_1<- - 1/(delta*log(1 - p))
  eta<- 1/delta
  # log-likelihood for discrete Pareto
  log_like<- function(N, delta, p){
   W1<- (1 - delta*(N - 1)*log(1 - p))^( - 1/delta)
   W2<- (1 - delta*N*log(1 - p))^( - 1/delta)
    return(sum(log( W1 - W2)+tol))
  }
  # Function to optimize to eta
  func_eta<-function(eta){
    ll<- eta*log(eta) - lgamma(eta)- eta - eta*log(mean(a)) + eta*mean(c)
    return( -ll)
  }
  # log-liklihood calcultion
  Devianceold<- 0
  Deviancenew <- log_like(N, delta, p)
  # Intitialize vectors for storing data
  Outi<- NULL; outd<- NULL; outp<-NULL; outD<- NULL; k = 1
  Outi[1]<- 0; outd[1]<- delta; outp[1]<- p; outD[1]<- Deviancenew
  while((abs(Deviancenew - Devianceold) > tol) & (k <= maxiter)){
    ### E step
    const<- 1/(gamma(eta)*((gamma_1 + N - 1)^(-eta) - (gamma_1 + N)^(-eta)))
    a<- const * gamma(eta + 1) * ((gamma_1 + N - 1)^(-(eta + 1)) - (gamma_1 + N)^(-(eta + 1)))
    t1<- digamma(eta) - log(gamma_1 + N - 1)
    t2<- digamma(eta) - log(gamma_1 + N)
    c<-  const*gamma(eta) * (t1*(gamma_1 + N - 1)^(-eta) - t2*(gamma_1 + N)^(-eta))
    #### M step
    constant<-mean(c)-log(mean(a))
    if(constant<0){
    eta <-stats:: nlm(f=func_eta,p=eta, ndigit = 12)$estimate
    }
    else{
      eta <- Inf
      }
    delta<- 1/eta
    p<- 1 - exp( - mean(a))
    gamma_1<- - 1/(delta*log(1 - p))
    Devianceold<-Deviancenew
    Deviancenew <- log_like(N, delta, p)
    # Output
    k<- k + 1
    Outi[k]<- k; outd[k]<- delta; outp[k]<- p; outD[k]<- Deviancenew
  }
  Output <- data.frame(Outi,outd,outp,outD)
  names(Output) <- c("iteration","delta","p","log-lik values")
  par<-data.frame(t(c(delta,p)))
  colnames(par)<-c("delta","p")
  result <- list(par=par, Deviance=Deviancenew, data=Output)
  return(result)
}

