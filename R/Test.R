#' Test of exponential versus Pareto type II (Lomax) distributions
#'
#' Performs \eqn{H_0:\alpha =0 \quad vs \quad H_0:\alpha =0}.  \eqn{\alpha =0} corresponds to exponential distribution and \eqn{\alpha >0} corresponds to Pareto type II distribution.
#'
#' expPareto_test is the function for testing exponential vs. Pareto type II distributions
#'
#' @param X data-vector.
#'
#' @return  vector of samples generate from BLG distribution.
#'
#' @examples
#' X<- Renext:: rlomax(100,scale = 0.200,shape = 0.5)
#' #X<-rexp(300, rate = 1/10)
#' test<- expPareto_test(X)
#' test
#'
#'@references  Arendarczyk, M. and Kozubowski T. J. and Panorska, A. k. (2018). A bivariate distribution with Lomax and geometric margins . Journal of the Korean Statistical Society, 47:405-422.
#' \url{https://doi.org/10.1016/j.jkss.2018.04.006}
#'
#' @export
expPareto_test <- function(X)
{
  test<- "Exponential vs. Pareto Type II test"
  data_name <- deparse(substitute(X))
  X <- X[!is.na(X)]
  n <- length(X)
  if (is.numeric(X))
    X <- as.double(X)
  else stop("argument 'X' must be numeric")
  ## statistic
  statistic <- sum(X^2)/(sum(X))^2
  #p-value calculation
  z <-  sqrt(n)*( ((n*statistic)/2) -1 )
  p_value<- stats:: pnorm(z,lower.tail = FALSE)
  ###Output
  result<-list(statistic=statistic, p.value=p_value, alternative="greater", data.name=data_name,
               method=test)
  class(result)<- "htest"
  return(result)
}
