#' EM algorithm function for the mixed bivariate distribution with mixture of exponential and geometric marginals (BMEG)
#'
#' This function computes the parameter estimates of BMEG distribution using the EM algorithm.
#'
#' Takes initial guess for the parameters and the algorithm will estmate the MLE.
#'
#' @param data  dataframe of bivariate random vector (X,N) BMEG distribution
#' @param beta scaled paramters vectors for the mixture components, which must be numeric greater than 0
#' @param q is vector of probabilities for belong to the components, which must between 0 and 1 and the sum equal to 1
#' @param K  number of mixture components
#' @param maxiter maximum number of iterations
#' @param tol tolerance value
#'
#' @return  list containg parameter estimates, Deviance and data frame of iteration
#'
#' @examples
#' data.df<-rbmixexpgeo(20, beta=c(1,2),p=0.6,q=c(0.4,0.6))
#' fit<-bmixexpgeo_em(data.df,K=2)
#' fit
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
   k = 1
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


#' The mixed bivariate distribution with mixture of exponential and geometric marginals (BMEG)
#'
#' Random sample generating function for the mixed bivariate distribution with mixture of exponential and geometric marginals (BMEG).
#'
#' rbexpgeo generates random sample from BEG distribution.
#'
#' @param n size of sample.
#' @param beta  vector of scale paramters which must be numeric greater than  0.
#' @param q vector of numeric parameters each between 0 and 1, and sum equal to 1.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of random samples generate from BMEG model.
#'
#' @examples
#' N<-rbmixexpgeo(20, beta=c(1,2),p=0.6,q=c(0.4,0.6))
#' N
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A Mixed bivariate distribution with mixture exponential and geometric marginals. Inprint.
#'
#' @export
rbmixexpgeo<- function(n,beta,p,q){
  N<- stats:: rgeom(n,p)+1
  u<- stats:: runif(n)
  q<-cumsum(q)
  indx<-findInterval(u,q)
  X<- stats:: rgamma(n,shape =N,rate = beta[c(indx)] )
  return(data.frame(X,N))
}


#' The mixed bivariate distribution with mixture of exponential and geometric marginals (BMEG)
#'
#' Density function for bivariate distribution with mixture of exponential and geometric marginals (BMEG).
#'
#' dbmixexpgeo is the density  function.
#'
#' @param data is bivariate vector  (X,N) vector representing observations from BMEG model.
#' @param beta  vector of scale paramters which must be numeric greater than  0.
#' @param q vector of numeric parameters each between 0 and 1, and sum equal to 1.
#' @param p numeric parameter between 0 and 1.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of densities.
#'
#' @examples
#' data<-rbmixexpgeo(20, beta=c(1,2),p=0.6,q=c(0.4,0.6))
#' den<-dbmixexpgeo(data, beta=c(1,2),p=0.6,q=c(0.4,0.6))
#' den
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A Mixed bivariate distribution with mixture exponential and geometric marginals. Inprint.
#'
#'
#' @export
dbmixexpgeo<- function(data,beta,p,q,log.p=FALSE){
  N<-data[,2]
  X<-data[,1]
  pdf_funct<-function(dat.df){
    pd<- q*(beta^dat.df[2])*(dat.df[1]^(dat.df[2]-1)) *
      exp(-beta*dat.df[1])*p*(1-p)^(dat.df[2]-1) /gamma(dat.df[2])
    return(pd)
  }
  M<-apply(data, 1, pdf_funct)
  if (log.p == FALSE){
    return(M)
  }else{
    M<- log(M)
    return(M)
  }
}


#' The mixed bivariate distribution with mixture of exponential and geometric marginals (BMEG)
#'
#' Distribution function for the mixed bivariate distribution with mixture of exponential and geometric marginals (BMEG).
#'
#' pbmixexpgeo is the distribution  function.
#'
#' @param data is bivariate vector  (X,N) vector representing observations from BMEG model.
#' @param beta  vector of scale paramters which must be numeric greater than  0.
#' @param q vector of numeric parameters each between 0 and 1, and sum equal to 1.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x, N \leq n]}, otherwise, \eqn{P[X > x, N > n]}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of distribution.
#'
#' @examples
#' data<-rbmixexpgeo(20, beta=c(1,2),p=0.6,q=c(0.4,0.6))
#' prob<-pbmixexpgeo(data, beta=c(1,2),p=0.6,q=c(0.4,0.6))
#' prob
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A Mixed bivariate distribution with mixture exponential and geometric marginals. Inprint.
#'
#' @export
pbmixexpgeo<- function(data,beta,p,q, lower.tail=TRUE,log.p=FALSE){
   ## cdf function
  cdf_funct<-function(dat.df){
  t <- seq(1,dat.df[2])
    cdf <- 0
    for (i in 1 :length(beta)){
      for(j in t){
        cdf<- cdf+ sum(pracma::gammainc(beta[i]*dat.df[1],dat.df[2]*j)[3] *q[i]*p*(1-p)^(j-1))
      }
    }

    return(cdf)
  }
  M<-apply(data, 1, cdf_funct)
  if (lower.tail==TRUE & log.p==FALSE){
    return(M)
  }
  else if (lower.tail==FALSE & log.p==TRUE){
    M<-log(1-M)
    return(M)
  }
  else if (lower.tail==TRUE & log.p==TRUE){
    M<-log(M)
    return(M)
  }
  else {
    return(1-M)
  }
}




