#' The Bivariate exponential geometric distrubution (BEG)
#'
#' BEG Parameter estimation
#'
#' bexpgeo_fit fits  BEG model to data.
#'
#' @param data is bivariate vector  (X,N) vector representing observations from BEG model
#' @param level confidence level espressed between 0 and 1 (Default is 0.95)
#'
#' @return  list of parameter estimates, confidence interval, deviance and covariance matrix
#'
#' @examples
#' Data.df<-rbexpgeo(n=100,beta = 10,p=0.45)
#' fit<-bexpgeo_fit(Data.df,level = 0.95)
#' fit
#'
#'@references  Kozubowski, T.J., & Panorska, A.K. (2005). A Mixed bivariate distribution with exponential and geometric marginals. Journal of Statistical Planning and Inference, 134, 501-520.
#' \url{https://doi.org/10.1016/j.jspi.2004.04.010}
#'
#' @export
bexpgeo_fit <- function(data,level=0.95) ## data has to be a vector (X,N)
{
  N<-data[,2] ## Get discrete data
  X<-data[,1] ## get continous data
  n<-length(N)
  qt<-(1-level)/2
  qz<- stats:: qnorm(qt)
  z<- abs(qz)
  b<-mean(N)/mean(X)
  p<-1/mean(N)
  sigmabb<-1/(b*b*p)
  sigmapp<-1/((1-p)*p*p)
  J= solve(matrix(c(sigmabb,0,0,sigmapp),byrow = 2, ncol = 2))
  lowerb<-b-z*sqrt(J[1,1]/n)
  lowerp<-p-z*sqrt(J[2,2]/n)
  upperb<-b+z*sqrt(J[1,1]/n)
  upperp<-p+z*sqrt(J[2,2]/n)
  log.like<-log(p)+ N*log(b)-lgamma(N)+(N-1)*log((1-p)*X)-b*X
  Deviance<- -2*sum(log.like)
  Output<-data.frame(matrix(c(b,p)),matrix(c(lowerb,lowerp)),matrix(c(upperb,upperp)))
  colnames(Output)<-c("estimate",paste(level*100,"%", " lower bound", sep=""),
                      paste(level*100,"%", " upper bound", sep=""))
  row.names(Output)<- c("beta","p")
  result <- list(Estimates=Output,Deviance=Deviance,Inverse.Fisher.Matrix=J)
  return(result)
}

