#' EM algorithm function for the mixed bivariate distribution with mixture of exponential and geometric marginals (BMEG)
#'
#' This function computes the parameter estimates of BMEG distribution using the EM algorithm.
#'
#' Takes initial guess for the parameters and the algorithm will estmate the MLE.
#'
#' @param data  dataframe of bivariate random vector (X,N) BMEG distribution
#' @param beta scaled paramters vectors for the mixture components, which must be numeric greater than 0
#' @param q is vector of probabilities for belong to the components, which must between 0 and 1 and the sum equal to 1
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
#'@references  Amponsah, C. K.,  Kozubowski, T. J. and Panorska (2019). A computational approach to estimation of discrete Pareto parameters. Inprint.
#'
#' @export
bmixexpgeo_em <- function(data, beta = NULL,q=NULL, K=2, maxiter = 500, tol = 1e-16){
  N <- data[,2]
  X<- data[,1]
  n <- length(N)
  # log-likelihood for BMEG distribution
  log_like<- function(data,beta,p,q){
    ## funct computes log of the sums of pdfs
    funct<-function(dat.df){
     pd<- q*(beta^dat.df[2])*(dat.df[1]^(dat.df[2]-1)) *
       exp(-beta*dat.df[1])*p*(1-p)^(dat.df[2]-1) /gamma(dat.df[2])
     return(log(sum(pd)))
    }
   ll<-apply(data, 1, funct)
    return(sum(ll))
  }
  #initialize beta and q
  if (beta==NULL && q==NULL){
    U<-runif(K)
    q<- U/sum(U)
  }
  #extimation of p
  p<-1/mean(N)
  # log-liklihood calcultion
  Devianceold<- 0
  Deviancenew <- log_like(data,beta,p,q)
  # Intitialize vectors for storing data
  Outi<- NULL; outb<- beta; outp<-NULL;outq<-NULL; outD<- NULL; k = 1
  Outi[1]<- 0; outb[1]<- beta; outp[1]<- p;outq<-q; outD[1]<- Deviancenew
  while((abs(Deviancenew - Devianceold) > tol) & (k <= maxiter)){
    ### E step
    ## funct computes log of the sums of pdfs
    tau_funct<-function(dat.df){
      pd<- q*(beta^dat.df[2])*(dat.df[1]^(dat.df[2]-1)) *
        exp(-beta*dat.df[1])*p*(1-p)^(dat.df[2]-1) /gamma(dat.df[2])
      return(pd/sum(pd))
    }
    tau<-apply(data, 1, tau_funct)
    #### M step
    q<-apply(tau, 1, mean)

    beta<-
    p<- p

    Devianceold<-Deviancenew
    Deviancenew <- log_like(data,beta,p,q)
    # Output
    k<- k + 1
  }
  result <- list(beta=par, q=q, p=p, log.like=Deviancenew, Iterations=k)
  return(result)
}

