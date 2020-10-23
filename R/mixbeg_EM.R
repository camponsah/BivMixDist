#' The finite mixture of bivariate distribution with mixture exponential and geometric (FMBEG) marginals.
#'
#' Random sample generating function for the FMBEG model.
#'
#' rfmbeg generates random sample from FMBEG distribution.
#'
#' @param n size of sample.
#' @param beta  vector of scale parameters which must be numeric greater than  0.
#' @param w vector of numeric weight parameters each between 0 and 1, and sum equal to 1.
#' @param p numeric parameter between 0 and 1.
#'
#' @return  vector of random samples generate from BMEG model.
#'
#' @examples
#' N <- rfmbeg(10, beta = c(1,2,10), p=0.5, w=c(0.3,0.2,0.5))
#' N
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A finite mixture of bivariate distribution with mixture exponential and geometric (FMBEG) marginals. Inprint.
#'
#' @export
rfmbeg<- function(n,beta,p,w){
  if (sum(w) != 1) stop("argument 'w' must sum to 1")
  if (length(w) != c(length(beta)*length(p))) stop("length of argument 'w' must be equal to the product of length of 'beta' and 'p'")
  pair_list <- expand.grid(beta, p)
  u<- stats:: runif(n)
  k <- 1:n
  data <- data.frame(k,u)
  q<-cumsum(w)
  rfun <- function(dat){
    indx <- findInterval(dat[2],q)+1
    para <- pair_list[indx,]
    N<- stats:: rgeom(1,prob = para$Var2)+1
    X<- stats:: rgamma(1,shape =N,rate = para$Var1)
    return(c(X,N))
  }
  M<-t(apply(data, 1, rfun))
  colnames(M) <- c("X", "N")
  return(data.frame(M))
  }


#' The finite mixture of bivariate distribution with mixture exponential and geometric (FMBEG) marginals.
#'
#' Density function for FMBEG model.
#'
#' dfmbeg is the density  function.
#'
#' @param data is bivariate vector  (X,N) vector representing observations from FMBEG model.
#' @param beta  vector of scale parameters which must be numeric greater than  0.
#' @param w vector of numeric weight parameters each between 0 and 1, and sum equal to 1.
#' @param p numeric parameters between 0 and 1.
#' @param log.p logical; if TRUE, densities are given as logarithmic values.
#'
#' @return  vector of densities.
#'
#' @examples
#' data <- rfmbeg(10, beta = c(1,2,10), p=0.5, w=c(0.3,0.2,0.5))
#' den <- dfmbeg(data, beta = c(1,2,10), p=0.5, w=c(0.3,0.2,0.5))
#' den
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A finite mixture of bivariate distribution with mixture exponential and geometric (FMBEG) marginals. Inprint.
#'
#'
#' @export
dfmbeg <- function(data,rate,p,w,log.p=FALSE){
  if (sum(w) != 1) stop("argument 'w' must sum to 1")
  if (length(w) != c(length(rate)*length(p))) stop("length of argument 'w' must be equal to the product of length of 'rate' and 'p'")
  N<-data[,2]
  X<-data[,1]
  expo <- dexp(X, rate = rate)
  dens <- function(dat.df){
    gam_pd <- stats:: dgamma(dat.df[1], shape = dat.df[2], rate = rate)
    geo_pd <- stats:: dgeom(dat.df[2], prob = p)/(1-p)
    pd <- sum(kronecker(gam_pd,geo_pd) *w)
    return(pd)
  }
  M <- apply(data, 1, dens)
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
#' @param beta  vector of scale parameters which must be numeric greater than  0.
#' @param w vector of numeric weight parameters each between 0 and 1, and sum equal to 1.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x, N \leq n]}, otherwise, \eqn{P[X > x, N > n]}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of distribution.
#'
#' @examples
#' data <- rfmbeg(10, beta = c(1,2,10), p=0.5, w=c(0.3,0.2,0.5))
#' prob <- pfmbeg(data, beta = c(1,2,10), p=0.5, w=c(0.3,0.2,0.5))
#' prob
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A Mixed bivariate distribution with mixture exponential and geometric marginals. Inprint.
#'
#' @export
pfmbeg <- function(data,beta,p,w, lower.tail=TRUE,log.p=FALSE){
  if (sum(w) != 1) stop("argument 'w' must sum to 1")
  if (length(w) != c(length(beta)*length(p))) stop("length of argument 'w' must be equal to the product of length of 'beta' and 'p'")
  pair_list <- data.frame(expand.grid(beta, p), w)
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
#' This function computes the parameter estimates of BMEG distribution using the EM algorithm.
#'
#' Takes initial guess for the parameters and the algorithm will estimate the MLE.
#'
#' @param data  dataframe of bivariate random vector (X,N) BMEG distribution.
#' @param beta scaled parameter vectors for the mixture components, which must be numeric greater than 0.
#' @param w is vector of probabilities for belong to the components, which must between 0 and 1 and the sum equal to 1.
#' @param m  number of gamma components
#' @param l  number of geometric components
#' @param maxiter maximum number of iterations.
#' @param tol tolerance value.
#' @param verb If TRUE, estimates are printed during each iteration.
#'
#' @return  list containing parameter estimates, Deviance and data frame of iteration
#'
#' @examples
#' data.df <- rfmbeg(10, beta = c(1,2,10), p=0.5, w=c(0.3,0.2,0.5))
#' fit <- fmbeg_em(data=data.df, m=3,l=1)
#' fit$par
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A Mixed bivariate distribution with mixture exponential and geometric marginals. Inprint.
#'
#' @export
fmbeg_em <- function(data, rate =NULL, p=NULL, w=NULL, m=1, l=2,
                     maxiter = 1000, tol = 1e-08,verb=FALSE){
  if( !is.numeric(m) | m <= 0) stop("Gamma components argument 'm' must numeric greater than '0'")
  if(!is.numeric(l) | l <= 0) stop("Geometric components argument 'l' must numeric greater than '0'")
  N <- data[,2]
  X <- data[,1]
  n <- length(N)
  if (is.null(w)) {
    beta.init <- fmbeg_gamma.init(X=data,rate = rate, w=NULL, m=m)
    p.init <- fmbeg_geo.init(X=data,p=p, w=NULL, l=l)
    rate <- beta.init$rate
    p <- p.init$p
    rate.w <- beta.init$w
    p.w <- p.init$w
    w <- kronecker(rate.w,p.w)/ sum(kronecker(rate.w,p.w))
  } else {
    rate.w <- w[1:m]/sum(w[1:m])
    p.w <- w[-c(1:m)]/sum(w[-c(1:m)])
  }
  pair_list <- data.frame(expand.grid(rate, p), w)
  colnames(pair_list) <- c("b","p","w")
 #dens<- function(d, df){
 #   pd<-NULL
 #   for(i in 1:nrow(d)){
 #    pd<- cbind(pd, dfmbeg(df, beta =d[i,]$b, p=d[i,]$p, w=d[i,]$w))
 #   }
 #   return(pd)
 # }
  ll.old <- sum(log(dfmbeg(data = data,rate = rate, p=p, w=w)))
  diff <- 1 + tol
  it <- 0
  ## Geometric mass function
  geo.dens <- function(par){
    d <- NULL
    for (i in 1:length(par)) {
      d <- cbind(d, par[i]*(1-par[i])^(N-1))
    }
    return(d)
  }
  ## Gamma density function
  gam.dens <- function(par){
    d <- NULL
    for (i in 1:length(par)) {
      d <- cbind(d, dgamma(X, shape = N, rate = par[i]))
    }
    return(d)
  }
  while(diff > tol && it < maxiter){
    ### Geometric membership probabilities
    z <- geo.dens(p)
    pi <- z/apply(z, 1, sum)
    p.w <- apply(pi, 2, mean)
    ### Gamma membership probabilities
    v <- gam.dens(rate)
    q <- v/apply(v, 1, sum)
    beta.w <- apply(pi, 2, mean)
    beta_hat<- apply(N *q, 2, mean)/apply(X *q, 2, mean)
    p_hat <- apply(pi, 2, mean)/apply(N *pi, 2, mean)
    w <- kronecker(rate.w,p.w)/ sum(kronecker(rate.w,p.w))
    pair_list <- data.frame(expand.grid(rate, p),w)
    colnames(pair_list) <- c("b","p","w")
    ll.new <- sum(log(dfmbeg(data = data,rate = rate, p=p, w=w)))
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
 para<- data.frame(expand.grid(rate, p),w)
 colnames(para) <- c("beta","p","w")
  result <- list(par=para, log.like=ll.new,
                 Iterations=it, posterior=q*pi, ft="fmbeg_em")

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
#' @param w is vector of probabilities for belong to the components, which must between 0 and 1 and the sum equal to 1
#' @param m  number of mixture components
#'
#' @return  list containing initial parameter estimates of beta and q
#'@references  Code adapted from Young et al. (2017). mixtools Package : Tools for Analyzing Finite Mixture Models, R CRAN.
#'
#' @export
fmbeg_gamma.init <- function(X, rate = NULL, w= NULL, m=2){
  n <- length(X[,1])
  if (is.null(w)) {
    u <- stats:: runif(m)
    w <- u/sum(u)
  }

  if(is.null(rate)){
  if(m==1){
    rate = mean(X[,2])/mean(X[,1])
  } else{
    x.sort<- X[order(X[,1]),]
    ind <- floor(n*cumsum(w))
    rate<- NULL
    t <- x.sort[(1:(ind[1]+1)),]
    rate[1] <- mean(t[,2])/mean(t[,1])
    for(j in 2:length(w)){
      t <- x.sort[(ind[j-1]:ind[j]),]
      rate[j] <- mean(t[,2])/mean(t[,1])
    }
  }
  ###
  }
 list(rate=rate, w=w)
}



#' parameter initialization for BMEG EM algorithm
#'
#' This function computes the parameter estimates of BMEG distribution using the EM algorithm.
#'
#' Takes initial guess for the parameters and the algorithm will estimate the MLE.
#'
#' @param N  discrete data from the marginal distribution of N in BMEG distribution
#' @param w is vector of probabilities for belong to the components, which must between 0 and 1 and the sum equal to 1
#' @param l  number of mixture components
#'
#' @return  list containing initial parameter estimates of beta and q
#'@references  Code adapted from Young et al. (2017). mixtools Package : Tools for Analyzing Finite Mixture Models, R CRAN.
#'
#' @export
fmbeg_geo.init <- function(X, p = NULL, w= NULL, l=1){
  N <- X[,2]
  n <- length(N)
  if (is.null(w)) {
    u <- stats:: runif(l)
    w <- u/sum(u)
  } #else k <- length(w)
  if(l==1){
    N.bar <- mean(N)
  } else{
    N.sort <- sort(N)
    ind <- floor(n*cumsum(w))
    N.part<-list()
    N.part[[1]] <- N.sort[1:(ind[1]+1)]
    for(j in 2:length(w)){
      N.part[[j]] <- N.sort[ind[j-1]:ind[j]]
    }
    N.bar <- sapply(N.part,mean)
  }
  if(is.null(p)){
    p <- 1/N.bar
    ind <- match(c(0,1),p)
    p <-  replace(p, ind, 0.5)
  }
  list(p=p, w=w)
}

