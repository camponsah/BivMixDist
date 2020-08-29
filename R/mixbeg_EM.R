#' The mixed bivariate distribution with mixture of exponential and geometric marginals (BMEG)
#'
#' Random sample generating function for the mixed bivariate distribution with mixture of exponential and geometric marginals (BMEG).
#'
#' rbexpgeo generates random sample from BEG distribution.
#'
#' @param n size of sample.
#' @param beta  vector of scale parameters which must be numeric greater than  0.
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
  indx<-findInterval(u,q)+1
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
#' @param beta  vector of scale parameters which must be numeric greater than  0.
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
  dens<-function(dat.df){
    pd<- q*(beta^dat.df[2])*(dat.df[1]^(dat.df[2]-1)) *
      exp(-beta*dat.df[1])*p*(1-p)^(dat.df[2]-1) /gamma(dat.df[2])
    return(pd)
  }
  M<-apply(data, 1, dens)
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
#' @param beta  vector of scale parameters which must be numeric greater than  0.
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


#' EM algorithm function for the mixed bivariate distribution with mixture of exponential and geometric marginals (BMEG)
#'
#' This function computes the parameter estimates of BMEG distribution using the EM algorithm.
#'
#' Takes initial guess for the parameters and the algorithm will estimate the MLE.
#'
#' @param data  dataframe of bivariate random vector (X,N) BMEG distribution.
#' @param beta scaled parameter vectors for the mixture components, which must be numeric greater than 0.
#' @param q is vector of probabilities for belong to the components, which must between 0 and 1 and the sum equal to 1.
#' @param k  number of mixture components
#' @param maxiter maximum number of iterations.
#' @param tol tolerance value.
#' @param verb If TRUE, estimates are printed during each iteration.
#'
#' @return  list containg parameter estimates, Deviance and data frame of iteration
#'
#' @examples
#' data.df<-rbmixexpgeo(20, beta=c(1,2),p=0.6,q=c(0.4,0.6))
#' fit<-bmixexpgeo_em(data.df,k=2)
#' fit$par
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A Mixed bivariate distribution with mixture exponential and geometric marginals. Inprint.
#'
#' @export
bmixexpgeo_em <- function(data, beta =NULL,q=NULL, k=2, maxiter = 1000,
                          tol = 1e-08,verb=FALSE){
  N <- data[,2]
  X<- data[,1]
  n <- length(N)
  fit <- bmixexpgeo.init(data,beta = beta, q=q, k=k)
  beta <- fit$beta
  q <- fit$q
  p<-1/mean(N)
  dens<- function(beta,q){
    pd<-NULL
    for(i in 1:k){
     pd<- cbind(pd, stats::dgamma(X,shape = N,rate = beta[i])* q[i]* p* (1-p)^(N-1))
    }
    return(pd)
  }
  ll.old <- sum(log(apply(dens(beta = beta,q=q),1,sum)))
  diff <- 1 + tol
  it <- 0
  while(diff > tol && it < maxiter){
    z <- dens(beta = beta,q=q)
    tau<- z/apply(z, 1, sum)
    q<-apply(tau, 2, mean)
    beta<- apply(N *tau, 2, mean)/apply(X *tau, 2, mean)
    ll.new <- sum(log(apply(dens(beta = beta,q=q),1,sum)))
    diff <- ll.new - ll.old
    ll.old <- ll.new
    it <- it +1
    if (verb) {
      cat("iteration =", it, " log.lik diff =", diff, " log.lik =",
          ll.new, "\n")
    }
  }
  if (it == maxiter) {
    cat("Warning! Convergence not achieved!", "\n")
  }
  beta <- rbind(beta)
  colnames(beta)=c(paste("comp", ".", 1:k, sep = ""))
  q <- rbind(q)
  colnames(q)=c(paste("comp", ".", 1:k, sep = ""))
  result <- list(data=data,par=list(beta=beta,q=q,p=p), log.like=ll.new,
                 Iterations=it, posterior=tau, ft="bmixexpgeo_em")

  class(result)<- "mixEM"
  return(result)
}



#' parameter initialization for BMEG EM algorithm
#'
#' This function computes the parameter estimates of BMEG distribution using the EM algorithm.
#'
#' Takes initial guess for the parameters and the algorithm will estimate the MLE.
#'
#' @param X  data-vector X from BMEG distribution
#' @param beta scale parameter vectors for the mixture components, which must be numeric greater than 0
#' @param q is vector of probabilities for belong to the components, which must between 0 and 1 and the sum equal to 1
#' @param k  number of mixture components
#'
#' @return  list containg initial parameter estimates of beta and q
#'@references  Code adapted from Young et al (2017). mixtools Package : Tools for Analyzing Finite Mixture Models, R CRAN.
#'
#' @export
bmixexpgeo.init <- function(X, beta = NULL, q = NULL, k=2){
  x <- X[,1]
  N <- X[,2]
  n <- length(x)

  if (is.null(q)) {
    U <- stats:: runif(k)
    q <- U/sum(U)
  } else k <- length(q)

  if(k==1){
    x.bar<- mean(x)
  } else{
    x.sort<- sort(x)
    ind=floor(n*cumsum(q))
    x.part<-list()
    x.part[[1]] <- x.sort[1:(ind[1]+1)]
    for(j in 2:k){
      x.part[[j]] <- x.sort[ind[j-1]:ind[j]]
    }
    x.bar <- sapply(x.part,mean)
  }
  if(is.null(beta)){
    beta <- mean(N)/x.bar
  }

  list( q=q, beta=beta, k=k)

}

