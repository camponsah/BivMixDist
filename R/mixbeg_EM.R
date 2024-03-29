#' The finite mixture of bivariate distribution with mixture exponential and geometric (FMBEG) marginals.
#'
#' Random sample generating function for the FMBEG model.
#'
#' rfmbeg generates random sample from FMBEG distribution.
#'
#' @param n size of sample.
#' @param beta  vector of parameters which must be numeric greater than  0.
#' @param p vector of numeric parameters between 0 and 1 of the finite mixture of geometric distribution.
#' @param q vector of gamma membership probabilities between 0 and 1, and sum equal to 1.
#' @param pi vector of geometric membership probabilities between 0 and 1, and sum equal to 1.
#'
#' @return  vector of random samples generate from BMEG model.
#'
#' @examples
#' data.df <- rfmbeg(10, beta = c(1,2,10), p=0.5, q=c(0.3,0.2,0.5), pi=1)
#' data.df
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A finite mixture of bivariate distribution with mixture exponential and geometric (FMBEG) marginals. Inprint.
#'
#' @export
rfmbeg<- function(n, beta, p, q, pi){
  if (sum(q) != 1) stop("argument 'q' must sum to 1")
  if (sum(pi) != 1) stop("argument 'pi' must sum to 1")
  if (c(length(q)*length(pi)) != c(length(beta)*length(p))) stop("Product of length of arguments 'q,pi' must be equal to the product of length of arguments 'beta, p'")
  pair_list <- NULL
  for(i in 1:length(q)){
    pair_list <- rbind(pair_list,expand.grid(q[i], pi))
  }
    if (length(p)==1){
      N <- stats:: rgeom(n,prob = p)+1
    } else{
      u <- stats:: runif(n)
      indx <- findInterval(u, cumsum(pi))+1
      indx <- data.frame(indx)
      fun <- function(par) {
        t <- stats:: rgeom(1, prob = p[par])+1
        return(t)
      }
      N <- apply(indx, 1, fun)
    }
  if(length(beta)==1){
    X <- stats:: rgamma(n,shape =N,rate = beta)
  } else{
    u <- stats:: runif(n)
    indx <- findInterval(u, cumsum(q))+1
    indx <- data.frame(indx,N)
    fun <- function(par) {
      t <- stats:: rgamma(1, shape =par[2], rate = beta[par[1]])
      return(t)
    }
    X <- apply(indx, 1, fun)
  }
  return(data.frame(X,N))
  }


#' The finite mixture of bivariate distribution with mixture exponential and geometric (FMBEG) marginals.
#'
#' Density function for FMBEG model.
#'
#' dfmbeg is the density  function.
#'
#' @param data is bivariate vector  (X,N) vector representing observations from FMBEG model.
#' @param beta  vector of parameters which must be numeric greater than  0.
#' @param p vector of numeric parameters between 0 and 1 of the finite mixture of geometric distribution.
#' @param q vector of gamma membership probabilities between 0 and 1, and sum equal to 1.
#' @param pi vector of geometric membership probabilities between 0 and 1, and sum equal to 1.
#' @param log.p logical; if TRUE, densities are given as logarithmic values.
#'
#' @return  vector of densities.
#'
#' @examples
#' data.df <- rfmbeg(10, beta = c(1,2,10), p=0.5, q=c(0.3,0.2,0.5), pi=1)
#' den <- dfmbeg(data=data.df, beta = c(1,2,10), p=0.5, q=c(0.3,0.2,0.5), pi=1)
#' den
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A finite mixture of bivariate distribution with mixture exponential and geometric (FMBEG) marginals. Inprint.
#'
#'
#' @export
dfmbeg <- function(data, beta, p, q, pi, log.p=FALSE){
  if (c(length(q)*length(pi)) != c(length(beta)*length(p))) stop("Product of length of arguments 'q,pi' must be equal to the product of length of arguments 'beta, p'")
  func_den <- function(df){
    den.gamm <- stats:: dgamma(x=df[1], shape = df[2], rate=beta) * q
    den.geo <- pi* p*(1-p)^(df[2]-1)
    den <- sum(kronecker(den.gamm, den.geo))
    return(den)
  }
  M <- apply(data, 1, func_den)
  if (log.p == FALSE){
    return(M)
  }else{
    M <- log(M)
    return(M)
  }
}


#' The finite mixture of bivariate distribution with mixture exponential and geometric (FMBEG) marginals.
#'
#' Distribution function for FMBEG.
#'
#' pfmbeg is the distribution  function.
#'
#' @param data is bivariate vector  (X,N) vector representing observations from BMEG model.
#' @param beta  vector of parameters which must be numeric greater than  0.
#' @param p vector of numeric parameters between 0 and 1 of the finite mixture of geometric distribution.
#' @param q vector of gamma membership probabilities between 0 and 1, and sum equal to 1.
#' @param pi vector of geometric membership probabilities between 0 and 1, and sum equal to 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x, N \leq n]}, otherwise, \eqn{P[X > x, N > n]}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of distribution.
#'
#' @examples
#' data.df <- rfmbeg(10, beta = c(1,2,10), p=0.5, q=c(0.3,0.2,0.5), pi=1)
#' prob <- pfmbeg(data=data.df, beta = c(1,2,10), p=0.5, q=c(0.3,0.2,0.5), pi=1)
#' prob
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A Mixed bivariate distribution with mixture exponential and geometric marginals. Inprint.
#'
#' @export
pfmbeg <- function(data, beta, p, q, pi, lower.tail=TRUE, log.p=FALSE){
  if(sum(q) != 1) stop("argument 'q' must sum to 1")
  if (sum(pi) != 1) stop("argument 'pi' must sum to 1")
  if (c(length(q)*length(pi)) != c(length(beta)*length(p))) stop("Product of length of arguments 'q,pi' must be equal to the product of length of arguments 'beta, p'")
  pair_list <- NULL
  for(i in 1:length(beta)){
    pair_list <- rbind(pair_list,expand.grid(beta[i], p))
  }
  w <- kronecker(q, pi)
  pair_list <- data.frame(pair_list, w)
  beg_cdf <- function(d){
    k<- 1:n
    df <- sum( stats:: pgamma(x, shape = k, rate = d[1])* d[2]*(1-d[2])^(k-1) ) *d[3]
    return(df)
  }
  M <- NULL
  for (i in nrow(data)) {
    x <- data[i,1]
    n <- data[i,2]
    M[i] <- sum(apply(pair_list , 1, beg_cdf))
  }
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
#' This function computes the parameter estimates of FMBEG distribution using the EM algorithm.
#'
#' Takes initial guess for the parameters and the algorithm will estimate the MLE.
#'
#' @param data  dataframe of bivariate random vector (X,N) BMEG distribution.
#' @param beta  numeric parameter vectors, which must be numeric greater than 0. Default value is NULL.
#' @param p vector of numeric parameters between 0 and 1. Default value is NULL.
#' @param q vector of gamma membership probabilities between 0 and 1, and sum equal to 1. Default value is NULL.
#' @param pi vector of geometric membership probabilities between 0 and 1, and sum equal to 1. Default value is NULL.
#' @param m  number of gamma components. Default value is 2.
#' @param l  number of geometric components. Default value is 1.
#' @param maxiter maximum number of iterations. Default value is 1000.
#' @param tol tolerance value. Default value is 1e-08
#' @param verb If TRUE, estimates are printed during each iteration.
#'
#' @return  list containing data frame parameter estimates, log-likelihood value and data frame of iteration
#'
#' @examples
#' data <- rfmbeg(10, beta = c(1,2,10), p=c(0.5,0.8), q=c(0.3,0.2,0.5), pi=c(0.6,0.4))
#' fit <- fmbeg_em(data=data, m=3,l=3)
#' fit$par
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A Mixed bivariate distribution with mixture exponential and geometric marginals. Inprint.
#'
#' @export
fmbeg_em <- function(data, beta =NULL, p=NULL, q=NULL, pi=NULL, m=2, l=1,
                     maxiter = 1000, tol = 1e-08,verb=FALSE){
  if( !is.numeric(m) | m <= 0) stop("Gamma components argument 'm' must numeric greater than '0'")
  if(!is.numeric(l) | l <= 0) stop("Geometric components argument 'l' must numeric greater than '0'")
  if (is.null(q) |is.null(beta)) {
    beta.init <- fmbeg_gamma.init(data, beta = beta, q=q, m=m)
    beta <- beta.init$beta
    q<- beta.init$q
  }
  if (is.null(pi) |is.null(p)) {
    p.init <- fmbeg_geo.init(data,p=p, pi=pi, l=l)
    p <- p.init$p
    pi <- p.init$pi
  }
    parameters <- function(beta,p,q,pi){
      pair_list <- NULL
      for(i in 1:length(beta)){
        pair_list <- rbind(pair_list,expand.grid(beta[i], p))
      }
      w <- kronecker(q,pi)
      pair_list <- data.frame(pair_list, w)
      colnames(pair_list) <- c("beta","p","weights")
      return(pair_list)
    }

  ll.old <- - sum(log(dfmbeg(data = data,beta = beta, p=p, q=q, pi=pi)))
  diff <- 1 + tol
  it <- 0
  ####
  while(diff > tol && it < maxiter){
    ### E & M steps
    fit_p <-  func_mgeo_em(data, p=p, pi=pi, l=l)
    fit_beta <- func_mgamma_em(data, beta=beta, q=q, m=m)
    ## M step
    pi <- fit_p$pi
    q <- fit_beta$q
    p <- fit_p$p
    beta <- fit_beta$beta
    ll.new <- - sum(log(dfmbeg(data = data,beta = beta, p=p, q=q, pi=pi)))
    diff <- abs(ll.new - ll.old)
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
 par <- parameters(beta,p,q,pi)
  result <- list(par=par, beta=beta, p=p, q=q, pi=pi, log.like=-ll.new,
                 Iterations=it, ft="fmbeg_em")

  class(result)<- "mixEM"
  return(result)
}



#' E and M steps for mxiture of geometric part of FMBEG
#'
#' This function computes the p and pi of FMBEG distribution.
#'
#' @param data  data-vector from FMBEG distribution.
#' @param p vector of numeric parameters between 0 and 1.
#' @param pi vector of geometric membership probabilities between 0 and 1, and sum equal to 1.
#' @param l  number of geometric components.
#'
#' @return  list containing  parameter estimates of p and pi
#'
#' @export
func_mgeo_em <- function(data, p, pi, l){
  dens_mgeo <- function(par){
    dens <- par[2]*par[1]*(1-par[1])^(data[,2]-1)
    return(dens)
  }
  par <- cbind(p,pi)
    ### E step
    z <- apply(par, 1, dens_mgeo)
    tau<- z/apply(z, 1, sum)
    ### M-step
    pi <- apply(tau, 2, mean)
    p <- apply(tau, 2, sum)/apply(data[,2]*tau, 2, sum)
  return(list(p=p, pi=pi))
}


#' E and M steps for mixture of gamma distributions part of FMBEG
#'
#' This function computes the beta and q of FMBEG distribution.
#'
#' @param data  data-vector from FMBEG distribution.
#' @param beta numeric parameters, which must be numeric greater than 0.
#' @param q vector of gamma membership probabilities between 0 and 1, and sum equal to 1
#' @param m  number of gamma components.
#'
#' @return  list containing  parameter estimates of p and pi
#'
#' @export
func_mgamma_em <- function(data, beta, q, m){
  dens_mgamma <- function(par){
    dens <- stats:: dgamma(x=data[,1], shape = data[,2], rate = par[1]) * par[2]
    return(dens)
  }
  par <- cbind(beta,q)
    ### E step
    z <- apply(par, 1, dens_mgamma)
    tau<- z/apply(z, 1, sum)
    ## M step
    q <- apply(tau, 2, mean)
    beta <- apply(data[,2]*tau, 2, sum)/apply(data[,1]*tau, 2, sum)
  return(list(beta=beta, q=q))
}



#' parameter initialization for BMEG EM algorithm
#'
#' This function computes the parameter estimates of BMEG distribution using the EM algorithm.
#'
#' Takes initial guess for the parameters and the algorithm will estimate the MLE.
#'
#' @param data  data-vector from FMBEG distribution
#' @param beta numeric parameters, which must be numeric greater than 0. Default value is NULL
#' @param q vector of gamma membership probabilities between 0 and 1, and sum equal to 1. Default value is NULL.
#' @param m  number of gamma components. Default value is 2.
#'
#' @return  list containing initial parameter estimates of beta and q
#'@references  Code adapted from Young et al. (2017). mixtools Package : Tools for Analyzing Finite Mixture Models, R CRAN.
#'
#' @export
fmbeg_gamma.init <- function(data, beta = NULL, q= NULL, m){
  n <- length(data[,1])
  if (is.null(q)) {
    u <- stats:: runif(m)
    q <- u/sum(u) #rep(1/m,m)
  }

  if(is.null(beta)){
  if(m==1){
    beta = mean(data[,2])/mean(data[,1])
  } else{
    x.sort<- data[order(data[,1]),]
    ind <- floor(n*cumsum(q))
    beta <- NULL
    t <- x.sort[(1:(ind[1]+1)),]
    beta[1] <- mean(t[,2])/mean(t[,1])
    for(j in 2:m){
      t <- x.sort[(ind[j-1]:ind[j]),]
      beta[j] <- mean(t[,2])/mean(t[,1])
    }
  }
  ###
  }
 return(list(beta=beta, q=q))
}



#' parameter initialization for BMEG EM algorithm
#'
#' This function computes the parameter estimates of BMEG distribution using the EM algorithm.
#'
#' Takes initial guess for the parameters and the algorithm will estimate the MLE.
#'
#' @param data  data from FMBEG distribution
#' @param p vector of numeric parameters between 0 and 1. Default value is NULL.
#' @param pi vector of geometric membership probabilities between 0 and 1, and sum equal to 1. Default value is NULL.
#' @param l  number of geometric components. Default value is 1.
#'
#' @return  list containing initial parameter estimates of p and membership probabilities
#'
#'@references  Code adapted from Young et al. (2017). mixtools Package : Tools for Analyzing Finite Mixture Models, R CRAN.
#'
#' @export
fmbeg_geo.init <- function(data, p = NULL, pi= NULL, l){
  N <- data[,2]
  n <- length(N)
  if (is.null(pi)) {
    u <- stats:: runif(l)
    pi <- u/sum(u)# rep(1/l, l)
  }
  if(l==1){
    N.bar <- mean(N)
  } else{
    N.sort <- sort(N)
    ind <- floor(n*cumsum(pi))
    N.part<-list()
    N.part[[1]] <- N.sort[1:(ind[1]+1)]
    for(j in 2:l){
      N.part[[j]] <- N.sort[ind[j-1]:ind[j]]
    }
    N.bar <- sapply(N.part,mean)
  }
  if(is.null(p)){
    p <- 1/N.bar
    ind <- match(c(0,1),p)
    p <-  replace(p, ind, 0.5)
  }
  list(p=p, pi=pi)
}

