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
#' N <- rfmbeg(10, beta = c(1,2,10), p=0.5, q=c(0.3,0.2,0.5), pi=1)
#' N
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
    }else{
      u <- stats:: runif(n)
      indx <- findInterval(u, cumsum(pi))+1
      indx <- data.frame(indx)
      fun <- function(par) {
        N <- stats:: rgeom(1, prob = p[par])+1
      }
      N <- apply(indx, 1, fun)
    }
  if(length(q)==1){
    X <- stats:: rgamma(n,shape =N,rate = beta)
  } else{
    u <- stats:: runif(n)
    indx <- findInterval(u, cumsum(q))+1
    indx <- data.frame(indx,N)
    fun <- function(par) {
      N <- stats:: rgamma(1, shape =par[2], rate = beta[par[1]])
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
#' data <- rfmbeg(10, beta = c(1,2,10), p=0.5, q=c(0.3,0.2,0.5), pi=1)
#' den <- dfmbeg(data, beta = c(1,2,10), p=0.5, q=c(0.3,0.2,0.5), pi=1)
#' den
#'
#'@references  Amponsah, C. K. and Kozubowski, T.J., and Panorska, A.K. (2020). A finite mixture of bivariate distribution with mixture exponential and geometric (FMBEG) marginals. Inprint.
#'
#'
#' @export
dfmbeg <- function(data, beta, p, q=1, pi=1, log.p=FALSE){
  if(sum(q) != 1) stop("argument 'q' must sum to 1")
  if (sum(pi) != 1) stop("argument 'pi' must sum to 1")
  if (c(length(q)*length(pi)) != c(length(beta)*length(p))) stop("Product of length of arguments 'q,pi' must be equal to the product of length of arguments 'beta, p'")
  dens <- function(dat.df){
    gam_pd <- q * stats:: dgamma(dat.df[1], shape = dat.df[2], rate = beta)
    geo_pd <- pi *  stats:: dgeom(dat.df[2], prob = p)/(1-p)
    pd <- prod(sum(gam_pd), sum(geo_pd))
    return(pd)
  }
  M <- apply(data.df, 1, dens)
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
#' data <- rfmbeg(10, rfmbeg(10, beta = c(1,2,10), p=0.5, q=c(0.3,0.2,0.5), pi=1)
#' prob <- pfmbeg(data, beta = c(1,2,10), p=0.5, q=c(0.3,0.2,0.5), pi=1)
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
fmbeg_em <- function(data, beta =NULL, p=NULL, q=NULL, pi=NULL, m=2, l=1,
                     maxiter = 1000, tol = 1e-08,verb=FALSE){
  if( !is.numeric(m) | m <= 0) stop("Gamma components argument 'm' must numeric greater than '0'")
  if(!is.numeric(l) | l <= 0) stop("Geometric components argument 'l' must numeric greater than '0'")
  N <- data[,2]
  X <- data[,1]
  n <- length(N)
  if (is.null(q) |is.null(beta)) {
    beta.init <- fmbeg_gamma.init(X=data, beta = beta, q=q, m=m)
    beta <- beta.init$beta
    q<- beta.init$q
  }
  if (is.null(pi) |is.null(p)) {
    p.init <- fmbeg_geo.init(X=data,p=p, pi=pi, l=l)
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
    ### parameters estimation function
    p_est <- function(data.df, beta, p, q, pi){
      if (length(p)==1){
        p_hat <- 1/mean(data.df[,2])
        pi_hat <- 1
      } else {
        p_hat <- NULL
        pi_hat <- NULL
        for (i in 1:length(p)) {
          p_mat <- NULL
          for (j in 1:length(beta)){
            p_mat <- cbind(p_mat, dfmbeg(data = data.df, beta = beta[j],p=p[i],q=q[j], pi=pi[i]))
          }
          p_mat <- p_mat/dfmbeg(data = data.df, beta = beta, p=p, q=q, pi=pi)
          p_hat[i] <- sum(apply(p_mat,2,sum))/sum(data.df[,2]*apply(p_mat,1,sum))
          pi_hat[i] <- sum(apply(p_mat,1,sum))
        }
        pi_hat <- pi_hat/sum(pi_hat)
        par<- cbind(p_hat, pi_hat)
        colnames(par) <- c("p", "pi")
        return(par)
      }
    }
    beta_est <- function(data.df, beta, p, q, pi){
        if (length(beta)==1){
          beta_hat <- mean(data.df[,2])/mean(data.df[,1])
          q_hat <- 1
        } else {
          beta_hat <- NULL
          q_hat <- NULL
          for (i in 1:length(beta)) {
            beta_hat <- NULL
            for (j in 1:length(p)){
              beta_hat <- cbind(beta_hat,dfmbeg(data = data, beta = beta[i],p=p[j], q=q[i],pi=pi[j]))
            }
            p_mat <- p_mat/dfmbeg(data = data.df, beta = beta, p=p, q=q, pi=pi)
            beta_hat[i]<- sum(data.df[,2]*apply(beta_hat,1,sum))/sum(data.df[,1]*apply(beta_hat,1,sum))
            q_hat[i] <- sum(apply(beta_mat,1,sum))
          }
          q_hat <- q_hat/sum(q_hat)
          par<- cbind(beta_hat, q_hat)
          colnames(par) <- c("beta", "q")
          return(par)
        }
    }
  ll.old <- sum(log(dfmbeg(data = data,beta = beta, p=p, q=q, pi=pi)))
  diff <- 1 + tol
  it <- 0
  ####
  while(diff > tol && it < maxiter){
    ### E and M steps
    fit_p <- p_est(data.df=data, beta=beta, p=p, q=q, pi=pi)
    fit_beta <- beta_est(data.df=data, beta=beta, p=p, q=q, pi=pi)
    beta <- fit_beta$beta
    q <- fit_beta$q
    p <- fit_p$p
    pi <- fit_p$pi
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
 par<- parameters(beta,p,q,pi)
  result <- list(par=para, beta=beta, p=p, q=q, pi=pi, log.like=ll.new,
                 Iterations=it, ft="fmbeg_em")

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
fmbeg_gamma.init <- function(X, beta = NULL, q= NULL, m=2){
  n <- length(X[,1])
  if (is.null(q)) {
    u <- stats:: runif(m)
    q <- u/sum(u)
  }

  if(is.null(beta)){
  if(m==1){
    beta = mean(X[,2])/mean(X[,1])
  } else{
    x.sort<- X[order(X[,1]),]
    ind <- floor(n*cumsum(w))
    beta<- NULL
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
#' @param N  discrete data from the marginal distribution of N in BMEG distribution
#' @param w is vector of probabilities for belong to the components, which must between 0 and 1 and the sum equal to 1
#' @param l  number of mixture components
#'
#' @return  list containing initial parameter estimates of beta and q
#'@references  Code adapted from Young et al. (2017). mixtools Package : Tools for Analyzing Finite Mixture Models, R CRAN.
#'
#' @export
fmbeg_geo.init <- function(X, p = NULL, pi= NULL, l=1){
  N <- X[,2]
  n <- length(N)
  if (is.null(pi)) {
    u <- stats:: runif(l)
    pi <- u/sum(u)
  } #else k <- length(w)
  if(l==1){
    N.bar <- mean(N)
  } else{
    N.sort <- sort(N)
    ind <- floor(n*cumsum(w))
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

