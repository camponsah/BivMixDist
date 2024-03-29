#' The multivariate geometric distribution (MGEO)
#'
#' Random sample generating function for multivariate geometric distribution with parameters $p_1,...,p_d$ in (0,1)$ and \eqn{-1<\theta < \dfrac{1-\min(p_1,...,p_d)}{1+\min(p_1,...,p_d)}}.
#'
#' rmgeo generates random sample from MGEO distribution.
#'
#' @param n size of sample.
#' @param prob vector of parameters between 0 and 1.
#' @param theta  numeric parameter which must be greater than -1 but less than (1-min(p))/(1+min(p)).
#'
#' @return  vector of random samples generate from MGEO model.
#'
#' @examples
#' N<-rmgeo(n=20, theta = 0.6, prob = c(0.2, 0.5))
#' N
#'
#'@references  Amponsah, C. K., & Kozubowski, T.J. (2022). A mixed multivariate distribution with bivariate exponential and geometric marginals. Inprint.
#' \url{https://doi.org/10.1016/j.jspi.2004.04.010}
#'
#' @export
rmgeo<- function(n, theta, prob){
  if ((n< 1) | !is.numeric(n)) {
    stop("Sample size n must be an integer whose value is at least 1!\n")
  }
  if (!is.numeric(prob) | sum(prob>=1) != 0 | sum(prob <= 0) != 0) {
    stop("Sample size n must be an integer whose value is at least 1!\n")
  }
  N <- stats:: rgeom(n, prob = prob[1]) +1
  if (length(prob)==1){
    return(N)
  }
  else {
    U <- stats:: runif(n)
    M <- NULL
    SF <- function(x) {
      u1 <- (1-prob[2])^x
      u2 <- (1+prob[2])^x
      a <- 1 - 2*(1/(1+prob[1]))^N[i]
      f <- u1 + theta* u1 *(1- 1/u2) *a - 1 + U[i]
      return(f)
    }
    for (i in 1 : n) {
      M[i] <- nleqslv::nleqslv(1, SF, jacobian=TRUE,control=list(btol=.01))$x
      M[i] <-  ceiling(M[i])
    }
    return(data.frame(N,M))
  }
}

#' The multivariate geometric distribution (MGEO)
#'
#' Probability mass function for multivariate geometric distribution with parameters $p_1,...,p_d$ in (0,1)$ and \eqn{-1<\theta < \dfrac{1-\min(p_1,...,p_d)}{1+\min(p_1,...,p_d)}}.
#'
#' dmgeo is the probability mass function.
#'
#' @param data is a data frame  MGEO model.
#' @param prob vector of parameters between 0 and 1.
#' @param theta  numeric parameter which must be greater than -1 but less than (1-min(p))/(1+min(p)).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of probabilities.
#'
#' @examples
#' data<-rmgeo(n=20, theta = 2, prob = c(0.6, 0.1))
#' den<-dmgeo(data, theta=2, prob = c(0.6, 0.1))
#' den
#'
#'@references  Kozubowski, T.J., & Panorska, A.K. (2005). A Mixed bivariate distribution with exponential and geometric marginals. Journal of Statistical Planning and Inference, 134, 501-520.
#' \url{https://doi.org/10.1016/j.jspi.2004.04.010}
#'
#' @export
dmgeo<- function(data, theta, prob, log.p=FALSE){
  w <- 1 - 2 * t(((1/(1+ prob)))^(t(data)))
  w <- 1 + theta * apply(w, 1, prod)
  t <- prod(prob)* t((1-prob)^(t(data) -1))
  t <- apply(t, 1, prod)
  M <- t * w
  if (log.p == FALSE){
    return(M)
  }else{
    M<- log(M)
    return(M)
  }
}


#' The bivariate exponential geometric distrubution (BEG)
#'
#' Distribution function for bivariate exponential geometric distribution with parameters \eqn{\beta > 0} and p in (0,1).
#'
#' pbexpgeo is the distribution  function.
#'
#' @param data is bivariate vector  (X,N) vector representing observations from BEG model.
#' @param beta  numeric parameter which must be numeric greater than 0.
#' @param p numeric parameter between 0 and 1.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x, N \leq n]}, otherwise, \eqn{P[X > x, N > n]}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return  vector of distribution.
#'
#' @examples
#' data<-rbexpgeo(20, beta=2,p=0.6)
#' den<-pbexpgeo(data, beta=2,p=0.6)
#' den
#'
#'@references  Kozubowski, T.J., & Panorska, A.K. (2005). A Mixed bivariate distribution with exponential and geometric marginals. Journal of Statistical Planning and Inference, 134, 501-520.
#' \url{https://doi.org/10.1016/j.jspi.2004.04.010}
#'
#' @export
pbexpgeo<- function(data,beta,p, lower.tail=TRUE,log.p=FALSE){
  N<-data[,2]
  X<-data[,1]
  S1<- stats:: ppois(N-1,lambda = ((1-p)*beta*X),lower.tail = TRUE)
  S2<- stats:: ppois(N-1,lambda = (beta*X),lower.tail = FALSE)
  M<-1-exp(-p*beta*X)*S1 -((1-p)^N)*S2
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



#' Fits the multivariate geometric distribution (MGEO)
#'
#' This function computes the parameter estimates of MGEO .
#'
#' mgeo_fit  fits MGEO model to data.
#'
#' @param data is data.frame of observations from MGEO model
#' @param tol is error tolerance level for checking convergence.
#'
#' @return  list of parameter estimates and deviance
#'
#' @examples
#' Data.df<- rmgeo(n = 100,theta = 10, prob= c(0.23, 0.7))
#' fit<-mgeo_fit(Data.df)
#' fit
#'
#'@references  Kozubowski, T.J., & Panorska, A.K. (2005). A Mixed bivariate distribution with exponential and geometric marginals. Journal of Statistical Planning and Inference, 134, 501-520.
#' \url{https://doi.org/10.1016/j.jspi.2004.04.010}
#'
#' @export
mgeo_fit <- function(data, tol = .Machine$double.eps^(1/2),
                     maxiter = 500) ## data has to be a dataframe (N,M)
{
  log.lik <- function(x){
    mn <- apply(data, 2, mean)
    ll1 <- log(x[1]) +log(x[2]) + (mn[1]-1)*log(1-x[1])
    ll1 <- ll1 + (mn[2]-1)*log(1-x[2])
    w1 <- 2* (1+x[1])^(- data[,1])
    w2 <- 2* (1+x[2])^(- data[,2])
    w <- (1-w1)*(1-w2)
    ll2 <- mean(log(1 + x[3] * w))
    return(-(ll1+ll2))
  }
  log.lik.theta <- function(x){
    w1 <- 2* (1 + p[1])^(- data[,1])
    w2 <- 2* (1 + p[2])^(- data[,2])
    w <- (1-w1)*(1-w2)
    w <- mean(log(1 + x * w))
    return(w)
  }
  p <- 1/apply(data, 2, mean)
  c <- (1 + min(p))/ (1-min(p))

  s <- seq(-1, c, (c + 1)/1000)
  ll.va <- apply(data.frame(s), 1, log.lik.theta)
  init_theta <- s[which(ll.va == max(ll.va))]
  if ( is.infinite(init_theta))
  {init_theta <- 0
  } else if ( length(init_theta) != 1)
  { init_theta <- 0
  }else {
    init_theta <- init_theta
  }
  par <- nloptr:: crs2lm(x0 = c(p, init_theta), fn =log.lik
               , lower= c(0,0,-1), upper =c(1,1,c)
               , pop.size = 100000
               , ranseed = NULL
               , xtol_rel = tol)$par
  p <- par[1:2]
  theta <- par[3]
  log.like<- sum(log(dmgeo(data = data, theta = theta
                           , prob=p)))
  Output<-t(data.frame(matrix(c(p, theta))))
  colnames(Output)<- c("p1", "p2", "theta")
  row.names(Output)<- NULL
  result <- list(par=Output,log.like=log.like)
   return(result)
}







