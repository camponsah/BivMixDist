


#' The multivariate mixed gamma-geometric distributions (MMGG)
#'
#' Density function for multivariate mixed gamma-geometric distribution with parameters \eqn{\alpha_1, \alpha_1 >0, \beta_1, \beta_1 > 0}, \eqn{p_1, P_2 \in (0,1)} and \eqn{-1\leq \theta \frac{1+\min(p_1,p_2)}{1+\min(p_1,p_2)}}.
#'
#' dmgammageo is the density function.
#'
#' @param data  bivariate vector  (X,N) observations from BGG model.
#' @param alpha  numeric vector of parameters which must be  greater than  0.
#' @param beta  numeric vector of parameters which must be  greater than  0.
#' @param prob numeric vector parameters between 0 and 1.
#' @param theta numeric  parameter between -1 and (1+min(p))/(1-min(p)).
#' @param log.p logical; if TRUE, densities are given as log(den).
#'
#' @return  vector of densities.
#'
#' @examples
#' data.df<-rmgammageo(20, alpha=c(0.5, 1.5), beta=c(20, 50), prob=c(0.6, 0.4), theta= -0.2)
#' den<-dmgammageo(data.df, alpha=c(0.5, 1.5), beta=c(20, 50), prob=c(0.6, 0.4), theta= -0.2)
#' den
#'
#'@references Amponsah, C. K. and Kozubowski, T.J., (2022). Inprint
#'
#' @export
dmgammageo<- function(data, alpha, beta, prob, theta, log.p=FALSE){
  den1 <- dbgammageo(data[, c(1,2)], alpha = alpha[1]
                     , beta = beta[1], p=prob[1])
  den2 <- dbgammageo(data[, c(3,4)], alpha = alpha[2]
                     , beta = beta[2], p=prob[2])
  den3 <- dmgeo(data = data[, c(2,4)], theta = theta, prob = prob)
  M <- den1 * den2 * den3
  if (log.p == FALSE){
    return(M)
  }else{
    M<- log(M)
    return(M)
  }
}



#' Fits the multivariate mixed gamma-geometric distributions (MMGG)
#'
#' This function computes the parameters of the model .
#'
#' mgammageo_fit  fits MMGG model to data.
#'
#' @param data is data.frame of observations, (X, N, Y, M)
#'
#' @return  list of parameter estimates and the log-likelihood value
#'
#' @examples
#' Data.df<- data.df<-rmgammageo(20, alpha=c(0.5, 1.5), beta=c(20, 50), prob=c(0.6, 0.4), theta= -0.2)
#' fit <- mgammageo_fit(Data.df)
#' fit
#'
#'@references Amponsah, C. K. and Kozubowski, T.J., (2022). Inprint
#'
#' @export
mgammageo_fit <- function(data) {
  N <- data[, 2]
  M <- data[, 4]
  X <- data[, 1]
  Y <- data[, 3]
  # Conditions to gurrantee MLEs of alphas
  cons1 <- log(mean(N)/mean(X)) + mean(N*log(X/N))/mean(N)
  cons2<- log(mean(M)/mean(Y)) + mean(M*log(Y/M))/mean(M)
  # functions to optimize
  if (cons1 >= 0) {
    stop("Warning! MLE of alpha1 do not exist!", "\n")
  } else if (cons2 >= 0) {
    stop("Warning! MLE of alpha2 do not exist!", "\n")
  }
  #
  log.lik.alpha1 <- function(par) {
    B<-par * mean(N)/mean(X)
    ll<- par * mean(N)*log(B+1e-16)-B*mean(X)+mean(par[1]*N*log(X)-lgamma(par[1]*N))
    return(-ll)
  }
  log.lik.alpha2 <- function(par) {
    C<-par[2]*mean(M)/mean(Y)
    ll<- par[2]*mean(M)*log(C+1e-16)-C*mean(Y)+mean(par[2]*M*log(Y)-lgamma(par[2]*M))
    return(-ll)
  }
  log.lik.p <- function(x){
    w <- 1 - 2 * t(((1+p))^(-t(data[, c(2, 4)])))
    w <- apply(w, 1, prod)
    l.like <- sum(log(1 + x*w))
    return(-l.like)
  }
  # estimate alphas and betas
  V1 <- stats:: var(X)
  V2 <- stats:: var(Y)
  alpha1 <- mean(N*((mean(X))^2)/V1)
  alpha2 <- mean(M*((mean(Y))^2)/V2)
  a1 <- stats:: nlm(log.lik.alpha1, p = alpha1)$estimate
  a2 <- stats:: nlm(log.lik.alpha2, p = alpha2)$estimate
  a <- c(a1, a2)
  b <- a * c( mean(N)/mean(X), mean(M)/mean(Y))
  # estimate p1, p2 and theta
  p <- 1/apply(data[, c(2, 4)], 2, mean)
  c <- (1 + min(p))/(1- min(p))
  theta <- stats::optimize(log.lik.p, interval = c(-1, c)
                           , tol = .Machine$double.eps^0.4)$minimum
  log.like<- sum(log(dmgammageo(data = data, alpha = a
                                , beta = b, prob=p, theta = theta)))
  par <- t(data.frame(c(a, b, p, theta)))
  colnames(par)<- c("alpha1", "alpha2", "beta1", "beta2", "p1", "p2", "theta")
  row.names(par) <- NULL
  result <- list(par=par, log.like=log.like)
  return(result)
}



